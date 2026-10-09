//! nlr-cli — command-line orchestration: wires config/seq/scan/assemble/output into a full pipeline.
//!
//! Built-in mot.txt / store.txt are embedded via `include_str!` and used by default when the
//! user does not pass `-x` / `-y`; explicit user paths take precedence.

use std::collections::HashMap;
use std::path::PathBuf;

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::{AnnotatorSignatureDefinition, DomainCategory, FlexibleSeedConfig};

/// Built-in mot.txt (PWM config), embedded at compile time for default distribution.
pub const EMBEDDED_MOT: &str = include_str!("../data/mot.txt");
/// Built-in store.txt (CDF config), embedded at compile time for default distribution.
pub const EMBEDDED_STORE: &str = include_str!("../data/store.txt");

/// Run configuration (derived from parsed CLI args).
pub struct RunConfig {
    pub input_fasta: PathBuf,
    /// mot.txt path; `None` -> use built-in embedded mot.
    pub mot_file: Option<PathBuf>,
    /// store.txt path; `None` -> use built-in embedded store.
    pub store_file: Option<PathBuf>,
    /// Fragment length (default 20000).
    pub fragment_length: usize,
    /// Overlap (default 2000).
    pub overlap: usize,
    /// Thread count (default auto-detected).
    pub threads: usize,
    /// Fragments per thread batch (default 1000, controls batch size).
    pub seqs_per_thread: usize,
    /// Checkpoint directory (saved after scan; reused to skip scan on rerun).
    pub checkpoint_dir: Option<PathBuf>,
    /// Assembly parameters.
    pub assemble: nlr_assemble::AssembleParams,
    /// Final motif-accept threshold (default 1e-5; higher values improve recall for divergent NLRs).
    pub motif_accept_p: f64,
    /// Motif prefilter threshold (default 1e-4).
    pub motif_prelim_p: f64,
    /// Domain-category declarations for motif IDs outside the built-in tables (e.g. "21=CC"),
    /// used by external motif libraries.
    pub motif_categories: Vec<String>,
    /// Extra seed combinations declared by an external library (e.g. "21,4"); repeatable.
    pub extra_seeds: Vec<String>,
    /// Extra signatures declared by an external library (e.g. "21,4").
    pub extra_signatures: Vec<String>,
    /// Flexible seeding config (category + rank based); `None` keeps exact-match seeding.
    pub flexible_seed: Option<FlexibleSeedConfig>,
}

impl RunConfig {
    /// Parse motif id lists such as "21,4" or "21,4,6"
    /// (comma-separated combinations).
    fn parse_id_lists(list: &[String]) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        for item in list {
            for combo in item.split(';') {
                let ids: Vec<u8> = combo
                    .split(',')
                    .filter_map(|x| x.trim().parse::<u8>().ok())
                    .collect();
                if ids.len() >= 2 {
                    out.push(ids);
                }
            }
        }
        out
    }

    /// Parse --motif-category "ID=CAT" entries (CAT in NBARC/LRR/TIR/CC/LINKER/NA).
    fn parse_motif_categories(list: &[String]) -> std::collections::HashMap<u8, DomainCategory> {
        let mut m = std::collections::HashMap::new();
        for item in list {
            for kv in item.split(',') {
                let kv = kv.trim();
                if kv.is_empty() {
                    continue;
                }
                if let Some((id, cat)) = kv.split_once('=') {
                    if let Ok(id) = id.trim().parse::<u8>() {
                        m.insert(id, DomainCategory::from_str(cat.trim().to_uppercase().as_str()));
                    }
                }
            }
        }
        m
    }

    pub fn new(input_fasta: PathBuf, mot_file: Option<PathBuf>, store_file: Option<PathBuf>) -> Self {
        RunConfig {
            input_fasta,
            mot_file,
            store_file,
            fragment_length: 20000,
            overlap: 2000,
            threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            seqs_per_thread: 1000,
            checkpoint_dir: None,
            assemble: nlr_assemble::AssembleParams::default(),
            motif_accept_p: 1e-5,
            motif_prelim_p: 1e-4,
            motif_categories: Vec::new(),
            extra_seeds: Vec::new(),
            extra_signatures: Vec::new(),
            flexible_seed: None,
        }
    }
}

/// Run result (for CLI output).
pub struct RunResult {
    pub nlrs: Vec<MotifList>,
    /// All motifs grouped by DNA sequence id (for `-m` / `-c` output).
    pub motifs_by_seq: HashMap<String, Vec<Motif>>,
    pub def: AnnotatorSignatureDefinition,
}

/// Core orchestration: input FASTA -> chop -> six-frame translation -> scan -> coord map -> assemble.
///
/// `progress` is an optional progress bar (incremented per batch).
pub fn run(config: &RunConfig) -> std::io::Result<RunResult> {
    run_with_progress(config, None)
}

/// Core orchestration with a progress bar.
pub fn run_with_progress(
    config: &RunConfig,
    progress: Option<&indicatif::ProgressBar>,
) -> std::io::Result<RunResult> {
    // Interrupt flag (SIGINT graceful shutdown).
    let interrupted = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
    let interrupted_flag = interrupted.clone();
    // Register SIGINT handler (once; duplicate registration errors are ignored).
    let _ = ctrlc::set_handler(move || {
        tracing::warn!("interrupt signal received, shutting down gracefully...");
        interrupted_flag.store(true, std::sync::atomic::Ordering::SeqCst);
    });

    // 1. Load config (Arc-shared across threads). User paths take precedence over built-in.
    let def_cfg = load_motif_definition(config)?;
    let mut signature_def = AnnotatorSignatureDefinition::new()
        .with_extra_categories(RunConfig::parse_motif_categories(&config.motif_categories))
        .with_extra_rules(RunConfig::parse_id_lists(&config.extra_seeds),
                          RunConfig::parse_id_lists(&config.extra_signatures));
    if let Some(cfg) = config.flexible_seed {
        signature_def = signature_def.with_flexible_seed(cfg);
    }
    let parser = nlr_scan::MotifParser::with_thresholds(
        def_cfg, config.motif_prelim_p, config.motif_accept_p);
    let parser = std::sync::Arc::new(parser);

    // 1.5 Checkpoint resume: if a checkpoint exists, assemble and return (skip scan).
    if let Some(dir) = &config.checkpoint_dir {
        if let Some(existing) = load_checkpoint(dir, signature_def.clone())? {
            return Ok(existing);
        }
    }

    // 2. Streaming chop + fragment-level parallel scan.
    //    Fragments are pulled from the chopper in fixed-size chunks and scanned in parallel
    //    (one rayon task per fragment), then the chunk is released — so peak memory is
    //    bounded by chunk_size rather than by genome size.
    let chopper = nlr_seq::SequenceChopper::from_file(
        &config.input_fasta,
        config.fragment_length,
        config.overlap,
    )?;
    let mut chopper = chopper;

    let num_threads = config.threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| std::io::Error::other(e.to_string()))?;

    let chunk_size = config.seqs_per_thread.max(1);
    let mut motifs_by_seq: HashMap<String, Vec<Motif>> = HashMap::new();
    let mut scanned: u64 = 0;

    loop {
        // Pull a chunk of fragments from the chopper (streaming).
        let mut chunk: Vec<nlr_seq::translate::BioSequence> = Vec::with_capacity(chunk_size);
        while chunk.len() < chunk_size {
            match chopper.next_sequence() {
                Some(f) => chunk.push(f),
                None => break,
            }
        }
        if chunk.is_empty() {
            break;
        }
        let chunk_len = chunk.len() as u64;

        // Fragment-level parallelism: each fragment is one rayon task (fine-grained load
        // balancing — a large contig no longer stalls the whole batch).
        let chunk_results: Vec<Vec<(String, Motif)>> = pool.install(|| {
            use rayon::prelude::*;
            chunk
                .par_iter()
                .filter_map(|fragment| {
                    if interrupted.load(std::sync::atomic::Ordering::SeqCst) {
                        return None;
                    }
                    Some(scan_one_fragment(fragment, &parser))
                })
                .collect()
        });

        // Merge into the global map immediately; the chunk is dropped after this iteration.
        for frags in chunk_results {
            for (seq_id, motif) in frags {
                motifs_by_seq.entry(seq_id).or_default().push(motif);
            }
        }

        scanned += chunk_len;
        if let Some(pb) = progress {
            pb.set_position(scanned);
        }
        if interrupted.load(std::sync::atomic::Ordering::SeqCst) {
            break;
        }
    }
    tracing::info!("chopping + scan done: {} fragments", scanned);

    if interrupted.load(std::sync::atomic::Ordering::SeqCst) {
        tracing::warn!("scan interrupted, outputting completed batch results");
    }

    // 3. Per-sequence sort + adjacent dedup (deterministic regardless of scan parallelism).
    for v in motifs_by_seq.values_mut() {
        v.sort();
        dedup_adjacent(v);
    }

    // 3.5. Signature 预过滤：按 (contig, strand, frame) 的局部 motif 簇判定，与切片边界无关。
    //      簇间距取 fragment_length（20 kb），与旧实现单个 fragment 的判定尺度对齐，
    //      但改由 motif 自身坐标定义簇，故不受切片相位影响。
    //      必须在排序 + 去重之后执行，因为簇的切分依赖排序结果。
    let cluster_gap = config.fragment_length as u64;
    for v in motifs_by_seq.values_mut() {
        filter_motifs_by_signature(v, &signature_def, cluster_gap);
    }

    // 4. Assemble NLRs per contig. Contigs are independent, so this is parallelized;
    //    results are keyed by sorted seq_id and flattened to keep output order deterministic
    //    (identical to the prior serial `for seq_id in sorted_seq_ids` loop).
    let mut seq_ids: Vec<&String> = motifs_by_seq.keys().collect();
    seq_ids.sort();
    let per_seq: Vec<(usize, Vec<MotifList>)> = pool.install(|| {
        use rayon::prelude::*;
        seq_ids
            .par_iter()
            .enumerate()
            .map(|(idx, seq_id)| {
                let motifs = motifs_by_seq.get(*seq_id).unwrap().clone();
                let assembled =
                    nlr_assemble::assemble(seq_id, motifs, &config.assemble, &signature_def);
                (idx, assembled)
            })
            .collect()
    });
    // Reorder by original sorted index → deterministic contig order.
    let mut per_seq = per_seq;
    per_seq.sort_by_key(|(idx, _)| *idx);
    let mut nlrs: Vec<MotifList> = Vec::new();
    for (_, mut assembled) in per_seq {
        nlrs.append(&mut assembled);
    }

    let result = RunResult {
        nlrs,
        motifs_by_seq,
        def: signature_def,
    };

    // 6. Save checkpoint (if configured).
    if let Some(dir) = &config.checkpoint_dir {
        let _ = save_checkpoint(dir, &result);
    }

    Ok(result)
}

/// Scan a single fragment: six-frame translate → find motifs → signature filter → DNA coord map.
/// Returns `(seq_id, motif)` pairs for this fragment.
fn scan_one_fragment(
    fragment: &nlr_seq::translate::BioSequence,
    parser: &nlr_scan::MotifParser,
) -> Vec<(String, Motif)> {
    let (id, offset) = parse_fragment_id(&fragment.identifier);
    let fragment_len = fragment.len() as u64; // actual length (last fragment may be < fragment_length)
    let protein_seqs = fragment.translate2protein();
    let mut out: Vec<(String, Motif)> = Vec::new();
    for pseq in &protein_seqs {
        let list = parser.find_motifs(&pseq.identifier, &pseq.sequence);
        if list.motifs.is_empty() {
            continue;
        }
        // 注意：signature 预过滤已从片段级上移到 (contig, strand, frame) 级，
        // 见 `filter_motifs_by_signature`。此处只负责扫描与坐标映射。
        let (frame, strand) = parse_frame(&pseq.identifier);
        for motif in list.motifs {
            let mut m = motif;
            // Map protein coordinates to genomic coordinates.
            // Use actual fragment length (not fixed fragment_length).
            m.set_dna(id.clone(), offset, fragment_len, frame, strand);
            out.push((id.clone(), m));
        }
    }
    out
}

/// Resolve mot/store sources (user path takes precedence over built-in) and load MotifDefinition.
fn load_motif_definition(config: &RunConfig) -> std::io::Result<nlr_config::MotifDefinition> {
    let (mot_str, store_str) = load_profile_text(&config.mot_file, &config.store_file)?;
    nlr_config::MotifDefinition::load_from_str(&mot_str, &store_str)
}

/// 读取 mot/store 文本：用户路径优先，否则回落到编译期内嵌的默认配置。
fn load_profile_text(
    mot_file: &Option<PathBuf>,
    store_file: &Option<PathBuf>,
) -> std::io::Result<(String, String)> {
    let mot = match mot_file {
        Some(p) => std::fs::read_to_string(p).map_err(|e| {
            std::io::Error::other(format!("cannot read mot file {}: {}", p.display(), e))
        })?,
        None => EMBEDDED_MOT.to_string(),
    };
    let store = match store_file {
        Some(p) => std::fs::read_to_string(p).map_err(|e| {
            std::io::Error::other(format!("cannot read store file {}: {}", p.display(), e))
        })?,
        None => EMBEDDED_STORE.to_string(),
    };
    Ok((mot, store))
}

/// Parse fragment id "{id}_{offset}" -> (id, offset).
fn parse_fragment_id(id: &str) -> (String, u64) {
    // offset is the number after the last '_'.
    if let Some(pos) = id.rfind('_') {
        let (head, tail) = id.split_at(pos);
        let offset = tail[1..].parse::<u64>().unwrap_or(0);
        (head.to_string(), offset)
    } else {
        (id.to_string(), 0)
    }
}

/// Parse "{id}_frame±n" -> (frame, strand).
fn parse_frame(id: &str) -> (u8, nlr_core::strand::Strand) {
    let frame_part = id.rsplit("_frame").next().unwrap_or("");
    let mut chars = frame_part.chars();
    let strand_char = chars.next().unwrap_or('+');
    let frame: u8 = chars.as_str().parse().unwrap_or(0);
    let strand = if strand_char == '-' {
        nlr_core::strand::Strand::Reverse
    } else {
        nlr_core::strand::Strand::Forward
    };
    (frame, strand)
}

/// Remove adjacent duplicates with the same motif id and DNA start.
///
/// 删除后不推进索引，因此连续三个及以上相同命中也能全部去掉（只保留第一条）；
/// 旧实现在删除后仍前进，遇到三个相同相邻 motif 会残留一条。
fn dedup_adjacent(motifs: &mut Vec<Motif>) {
    let mut i = 0;
    while i < motifs.len() {
        if i + 1 < motifs.len()
            && motifs[i].id == motifs[i + 1].id
            && motifs[i].dna_start == motifs[i + 1].dna_start
        {
            motifs.remove(i + 1);
        } else {
            i += 1;
        }
    }
}

/// 局部 motif 簇级的 signature 预过滤（与切片边界无关）。
///
/// 背景：旧实现把预过滤放在单个 fragment（约 20 kb 窗口、单帧）上，匹配不到 signature 就把该
/// fragment 的全部命中丢掉。signature 要求 motif 邻接，因此切片边界会切断邻接关系，结果随
/// overlap / fragment_length 变化（实测同一区间在不同边界下 motif 数不同，个别位点丢失 LRR）。
///
/// 现改为按 (contig, strand, frame) 分组、再按 motif 自身坐标切成局部簇：
/// - 分组键只依赖 motif 的链与读码框，簇边界只依赖相邻 motif 的间距，与切片边界无关；
/// - `motifs` 已按「链优先 + 坐标」排序，组内顺序即该读码框的翻译顺序（正向升序、反向降序）；
/// - 相邻 motif 间距超过 `cluster_gap` 即断开新簇，等价于把「一个局部候选位点区域」当作判定单元，
///   与旧实现的 20 kb 窗口尺度接近但相位固定，故能保留过滤强度又消除边界依赖。
fn filter_motifs_by_signature(
    motifs: &mut Vec<Motif>,
    def: &AnnotatorSignatureDefinition,
    cluster_gap: u64,
) {
    if motifs.is_empty() {
        return;
    }
    // 6 个桶：正向 0/1/2 帧 + 反向 0/1/2 帧；桶内索引保持排序后的翻译顺序。
    let mut buckets: [Vec<usize>; 6] = Default::default();
    for (i, m) in motifs.iter().enumerate() {
        let base = if m.strand == nlr_core::strand::Strand::Forward { 0 } else { 3 };
        buckets[base + (m.frame as usize % 3)].push(i);
    }

    let mut keep = vec![false; motifs.len()];
    for bucket in buckets.iter() {
        if bucket.is_empty() {
            continue;
        }
        // 沿翻译顺序把该帧的 motif 切成局部簇（间距 > cluster_gap 断开）。
        let mut start = 0usize;
        while start < bucket.len() {
            let mut end = start + 1;
            while end < bucket.len() {
                let prev = &motifs[bucket[end - 1]];
                let cur = &motifs[bucket[end]];
                if cur.dna_start.abs_diff(prev.dna_start) > cluster_gap {
                    break;
                }
                end += 1;
            }
            // 簇内命中任一 signature（或满足柔性播种）则保留整簇。
            let ids: Vec<u8> = bucket[start..end].iter().map(|&i| motifs[i].id).collect();
            let matched = match def.flexible_seed() {
                Some(cfg) => def.has_signature_flexible(&ids, &cfg),
                None => def.has_signature(&ids),
            };
            if matched {
                for &i in &bucket[start..end] {
                    keep[i] = true;
                }
            }
            start = end;
        }
    }

    let mut idx = 0usize;
    motifs.retain(|_| {
        let k = keep[idx];
        idx += 1;
        k
    });
}

/// Flatten all motifs into a list (for `-m` / `-c` output).
pub fn all_motifs(result: &RunResult) -> Vec<Motif> {
    let mut out = Vec::new();
    for v in result.motifs_by_seq.values() {
        for m in v {
            out.push(m.clone());
        }
    }
    out
}

/// Save checkpoint: write motif results in `-c` TSV format into the checkpoint directory.
pub fn save_checkpoint(dir: &std::path::Path, result: &RunResult) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)?;
    let path = dir.join("motifs.tsv");
    let mut f = std::fs::File::create(&path)?;
    nlr_output::export_motifs(&mut f, &all_motifs(result))?;
    tracing::info!("checkpoint saved: {}", path.display());
    Ok(())
}

/// Load checkpoint: read motif results from the checkpoint directory.
pub fn load_checkpoint(
    dir: &std::path::Path,
    def: AnnotatorSignatureDefinition,
) -> std::io::Result<Option<RunResult>> {
    let path = dir.join("motifs.tsv");
    if !path.exists() {
        return Ok(None);
    }
    let text = std::fs::read_to_string(&path)?;
    let motifs = nlr_output::import_motifs(&text);
    let mut motifs_by_seq: HashMap<String, Vec<Motif>> = HashMap::new();
    for m in motifs {
        let seq = m.dna_sequence_id.clone().unwrap_or_default();
        motifs_by_seq.entry(seq).or_default().push(m);
    }
    for v in motifs_by_seq.values_mut() {
        v.sort();
        dedup_adjacent(v);
    }
    // Assemble (reuse existing motif results, skip scan).
    let mut nlrs: Vec<MotifList> = Vec::new();
    let mut seq_ids: Vec<&String> = motifs_by_seq.keys().collect();
    seq_ids.sort();
    for seq_id in seq_ids {
        let ms = motifs_by_seq.get(seq_id).unwrap().clone();
        let mut assembled = nlr_assemble::assemble(seq_id, ms, &nlr_assemble::AssembleParams::default(), &def);
        nlrs.append(&mut assembled);
    }
    tracing::info!("checkpoint loaded: {} sequences", motifs_by_seq.len());
    Ok(Some(RunResult {
        nlrs,
        motifs_by_seq,
        def,
    }))
}

// =====================================================================================
// 蛋白模式（--protein）：HMMER 结构域扫描 + PWM LRR 判定，输出独立结果集（无基因组坐标）
// =====================================================================================

/// 蛋白模式运行配置。
pub struct ProteinConfig {
    pub input_fasta: PathBuf,
    /// 用户追加/覆盖的 HMM 文件（可含多个模型）。
    pub hmm_files: Vec<PathBuf>,
    pub mot_file: Option<PathBuf>,
    pub store_file: Option<PathBuf>,
    pub threads: usize,
    /// 每批处理的蛋白条数（控制内存与进度粒度）。
    pub batch: usize,
    /// LRR PWM 的最终接受阈值。
    pub motif_accept_p: f64,
    /// LRR PWM 的预筛阈值。
    pub motif_prelim_p: f64,
}

impl ProteinConfig {
    pub fn new(input_fasta: PathBuf) -> Self {
        ProteinConfig {
            input_fasta,
            hmm_files: Vec::new(),
            mot_file: None,
            store_file: None,
            threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            batch: 256,
            motif_accept_p: 1e-5,
            motif_prelim_p: 1e-4,
        }
    }
}

/// 单条蛋白的判定结果。
pub struct ProteinRecord {
    pub id: String,
    pub length: usize,
    pub arch: nlr_domain::NlrArchitecture,
    /// 全部被 HMMER 报告的结构域命中（含未通过 GA 的，写入 domains.tsv）。
    pub hits: Vec<nlr_domain::DomainHit>,
    /// PWM 的 LRR 命中：(motif id, start, end, pvalue)，坐标均为蛋白上的 1-based 闭区间。
    pub lrr_sites: Vec<(u8, usize, usize, f64)>,
}

impl ProteinRecord {
    /// 通过 GA 的 HMM 命中，按位置排序。
    pub fn passed_domains(&self) -> Vec<&nlr_domain::DomainHit> {
        self.hits.iter().filter(|h| h.passed).collect()
    }
}

/// 蛋白模式整体结果。
pub struct ProteinResult {
    /// 全部输入蛋白的判定，按 id 排序（并行扫描后统一排序，保证输出可复现）。
    pub records: Vec<ProteinRecord>,
    /// 内置 + 用户模型的展示串，例如 `PF00931(NB-ARC)`。
    pub model_labels: Vec<String>,
}

/// 蛋白模式运行元信息（用于写摘要）。
pub struct ProteinRunMeta {
    pub input_name: String,
    pub threads: usize,
    pub elapsed: std::time::Duration,
    pub lrr_threshold: f64,
    pub hmm_files: Vec<String>,
}

/// 蛋白模式主流程：读蛋白 FASTA → 并行结构域扫描 → 汇总。
pub fn run_protein(
    config: &ProteinConfig,
    progress: Option<&indicatif::ProgressBar>,
) -> std::io::Result<ProteinResult> {
    // 1) 模型库：内置 4 个精选 Pfam 模型 + 用户 --hmm 覆盖/追加。
    let mut library = nlr_domain::HmmLibrary::from_embedded();
    for path in &config.hmm_files {
        let added = library.add_file(path)?;
        tracing::info!("loaded {added} domain model(s) from {}", path.display());
    }
    let model_labels: Vec<String> = library
        .models()
        .iter()
        .map(|m| format!("{}({})", m.acc, m.name))
        .collect();

    // 2) LRR 用现有 PWM（motif 9/11），与 DNA 模式共用同一份表与阈值。
    let (mot, store) = load_profile_text(&config.mot_file, &config.store_file)?;
    let definition = nlr_config::MotifDefinition::load_from_str(&mot, &store)?;
    let parser = std::sync::Arc::new(nlr_scan::MotifParser::with_thresholds(
        definition,
        config.motif_prelim_p,
        config.motif_accept_p,
    ));
    let signature_def = nlr_core::signature_def::AnnotatorSignatureDefinition::new();

    // 3) 流式读蛋白 + 批量并行扫描（扫描器在工作线程内构造，OProfile 不是 Send）。
    let mut reader = nlr_seq::fasta::FastaReader::from_file(&config.input_fasta)?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads.max(1))
        .build()
        .map_err(|e| std::io::Error::other(e.to_string()))?;

    let batch_size = config.batch.max(1);
    let mut records: Vec<ProteinRecord> = Vec::new();
    let mut scanned = 0u64;
    loop {
        let mut batch: Vec<nlr_seq::translate::BioSequence> = Vec::with_capacity(batch_size);
        while batch.len() < batch_size {
            match reader.read_entry() {
                Some(seq) => batch.push(seq),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }

        let batch_records: Vec<ProteinRecord> = pool.install(|| {
            use rayon::prelude::*;
            batch
                .par_iter()
                .map_init(
                    || library.scanner(),
                    |scanner, seq| scan_one_protein(scanner, &parser, &signature_def, seq),
                )
                .collect()
        });
        scanned += batch_records.len() as u64;
        records.extend(batch_records);
        if let Some(pb) = progress {
            pb.set_position(scanned);
        }
    }

    records.sort_by(|a, b| a.id.cmp(&b.id));
    tracing::info!("protein scan complete: {scanned} sequences");
    Ok(ProteinResult {
        records,
        model_labels,
    })
}

/// 扫描单条蛋白：HMM 结构域 + PWM LRR → 架构判定。
fn scan_one_protein(
    scanner: &mut nlr_domain::Scanner,
    parser: &nlr_scan::MotifParser,
    signature_def: &nlr_core::signature_def::AnnotatorSignatureDefinition,
    seq: &nlr_seq::translate::BioSequence,
) -> ProteinRecord {
    let hits = scanner.scan_protein(&seq.identifier, &seq.sequence);

    // PWM 只取 LRR 类 motif（默认 motif_9 / motif_11），find_motifs 内部已按 accept 阈值过滤。
    let list = parser.find_motifs(&seq.identifier, &seq.sequence);
    let mut lrr_sites: Vec<(u8, usize, usize, f64)> = list
        .motifs
        .iter()
        .filter(|m| signature_def.is_lrr(m.id))
        .map(|m| {
            let start = m.position as usize;
            let end = start + m.protein_sequence.len().saturating_sub(1);
            (m.id, start, end, m.pvalue)
        })
        .collect();
    lrr_sites.sort_by_key(|(_, s, e, _)| (*s, *e));

    let arch = nlr_domain::classify(&hits, lrr_sites.len());
    ProteinRecord {
        id: seq.identifier.clone(),
        length: seq.sequence.len(),
        arch,
        hits,
        lrr_sites,
    }
}

/// 写出蛋白模式结果集（domains.tsv / nlr.tsv / summary.txt / plots）。
pub fn write_protein_outputs(
    output_dir: &std::path::Path,
    prefix: &str,
    result: &ProteinResult,
    meta: &ProteinRunMeta,
) -> std::io::Result<Vec<PathBuf>> {
    use std::collections::BTreeMap;
    use std::io::Write as _;

    std::fs::create_dir_all(output_dir)?;
    let mut written: Vec<PathBuf> = Vec::new();

    // 结构域明细：HMM 命中 + PWM LRR 命中。
    let p = output_dir.join(format!("{prefix}.domains.tsv"));
    let mut f = std::fs::File::create(&p)?;
    writeln!(
        f,
        "#protein_id\tsource\tmodel\tmodel_name\tcategory\tstart\tend\tbitscore\tpvalue\tpassed"
    )?;
    for r in &result.records {
        for h in &r.hits {
            writeln!(
                f,
                "{}\thmm\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
                r.id,
                h.model_acc,
                h.model_name,
                h.category.as_str(),
                h.start,
                h.end,
                h.bitscore,
                nlr_core::motif::format_double_java(h.pvalue),
                h.passed
            )?;
        }
        for (id, start, end, pvalue) in &r.lrr_sites {
            writeln!(
                f,
                "{}\tpwm\tmotif_{}\tLRR\tLRR\t{}\t{}\t.\t{}\ttrue",
                r.id,
                id,
                start,
                end,
                nlr_core::motif::format_double_java(*pvalue)
            )?;
        }
    }
    written.push(p);

    // NLR 候选汇总。
    let p = output_dir.join(format!("{prefix}.nlr.tsv"));
    let mut f = std::fs::File::create(&p)?;
    writeln!(
        f,
        "#protein_id\tlength\tclass\tnterm\thas_nbarc\thas_tir\thas_lrr\tcomplete\tn_domains\tdomains"
    )?;
    let mut type_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut model_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut nlr_total = 0usize;
    let mut complete_total = 0usize;
    for r in &result.records {
        let passed = r.passed_domains();
        for h in &passed {
            *model_counts.entry(h.model_name.clone()).or_default() += 1;
        }
        if !r.lrr_sites.is_empty() {
            *model_counts.entry("LRR(PWM)".to_string()).or_default() += 1;
        }
        if !r.arch.is_nlr {
            continue;
        }
        nlr_total += 1;
        if r.arch.complete {
            complete_total += 1;
        }
        *type_counts.entry(r.arch.class.clone()).or_default() += 1;

        // 结构域摘要：HMM 通过 GA 的命中 + LRR 位点，按起点排序。
        let mut items: Vec<(usize, String)> = passed
            .iter()
            .map(|h| (h.start, format!("{}:{}-{}", h.model_name, h.start, h.end)))
            .collect();
        items.extend(
            r.lrr_sites
                .iter()
                .map(|(id, s, e, _)| (*s, format!("motif_{id}:{s}-{e}"))),
        );
        items.sort();
        let domains = items
            .iter()
            .map(|(_, s)| s.clone())
            .collect::<Vec<_>>()
            .join(";");
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.id,
            r.length,
            r.arch.class,
            r.arch.nterm,
            r.arch.has_nbarc,
            r.arch.has_tir,
            r.arch.has_lrr,
            r.arch.complete,
            items.len(),
            domains
        )?;
    }
    written.push(p);

    // 统计图（复用 nlr-plot 的通用计数柱状图）。
    let plot_dir = output_dir.join(format!("{prefix}.plots"));
    std::fs::create_dir_all(&plot_dir)?;
    let p1 = plot_dir.join("01-nlr-types.png");
    if let Err(e) = nlr_plot::plot_counts(&p1, "NLR types (protein mode)", "protein count", &type_counts)
    {
        tracing::warn!("NLR type plot generation failed: {e}");
    }
    written.push(p1);
    let p2 = plot_dir.join("02-domain-counts.png");
    if let Err(e) = nlr_plot::plot_counts(&p2, "Domain hits (protein mode)", "hit count", &model_counts)
    {
        tracing::warn!("domain count plot generation failed: {e}");
    }
    written.push(p2);

    // 运行摘要。
    use std::fmt::Write as _;
    let mut s = String::new();
    let _ = writeln!(s, "# FastNLR protein-mode run summary");
    let _ = writeln!(s, "# input file\t{}", meta.input_name);
    let _ = writeln!(s, "# threads\t{}", meta.threads);
    let _ = writeln!(s, "# elapsed\t{:.2}s", meta.elapsed.as_secs_f64());
    let _ = writeln!(s, "# domain models\t{}", result.model_labels.join(", "));
    if !meta.hmm_files.is_empty() {
        let _ = writeln!(s, "# user HMM files\t{}", meta.hmm_files.join(", "));
    }
    let _ = writeln!(
        s,
        "# LRR model\tPWM motif_9/motif_11 (p < {})",
        meta.lrr_threshold
    );
    let _ = writeln!(s, "# proteins scanned\t{}", result.records.len());
    let _ = writeln!(s, "# NLR candidates\t{}", nlr_total);
    let _ = writeln!(s, "# complete (NB-ARC + LRR)\t{}", complete_total);
    let _ = writeln!(s, "# NLR types (count desc):");
    let mut types: Vec<(&String, &usize)> = type_counts.iter().collect();
    types.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    for (class, count) in types {
        let _ = writeln!(s, "# {}\t{}", count, class);
    }
    let _ = writeln!(s, "# output files:");
    for path in &written {
        let _ = writeln!(s, "# {}", path.display());
    }
    let p = output_dir.join(format!("{prefix}.summary.txt"));
    std::fs::write(&p, &s)?;
    written.push(p);

    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nlr_core::strand::Strand;

    fn dna_motif(id: u8, start: u64) -> Motif {
        let mut m = Motif::new_protein(id, "chr1".to_string(), 1, "AAA".to_string(), 1e-6);
        m.set_dna("chr1".to_string(), start, 20_000, 0, Strand::Forward);
        m
    }

    #[test]
    fn dedup_adjacent_removes_all_consecutive_duplicates() {
        // 连续三个及以上相同命中（同 id 同坐标）应全部收敛为一条。
        let mut motifs = vec![
            dna_motif(6, 100),
            dna_motif(6, 100),
            dna_motif(6, 100),
            dna_motif(4, 200),
            dna_motif(4, 200),
            dna_motif(6, 300),
        ];
        dedup_adjacent(&mut motifs);
        let keys: Vec<(u8, u64)> = motifs.iter().map(|m| (m.id, m.dna_start)).collect();
        assert_eq!(keys, vec![(6, 100), (4, 200), (6, 300)]);
    }
}
