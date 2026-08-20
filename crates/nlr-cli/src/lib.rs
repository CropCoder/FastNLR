//! nlr-cli — command-line orchestration: wires config/seq/scan/assemble/output into a full pipeline.
//!
//! Built-in mot.txt / store.txt are embedded via `include_str!` and used by default when the
//! user does not pass `-x` / `-y`; explicit user paths take precedence.

use std::collections::HashMap;
use std::path::PathBuf;

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;

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
    /// Overlap (default 5000).
    pub overlap: usize,
    /// Thread count (default auto-detected).
    pub threads: usize,
    /// Fragments per thread batch (default 1000, controls batch size).
    pub seqs_per_thread: usize,
    /// Checkpoint directory (saved after scan; reused to skip scan on rerun).
    pub checkpoint_dir: Option<PathBuf>,
    /// Assembly parameters.
    pub assemble: nlr_assemble::AssembleParams,
}

impl RunConfig {
    pub fn new(input_fasta: PathBuf, mot_file: Option<PathBuf>, store_file: Option<PathBuf>) -> Self {
        RunConfig {
            input_fasta,
            mot_file,
            store_file,
            fragment_length: 20000,
            overlap: 5000,
            threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            seqs_per_thread: 1000,
            checkpoint_dir: None,
            assemble: nlr_assemble::AssembleParams::default(),
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
    let signature_def = AnnotatorSignatureDefinition::new();
    let parser = nlr_scan::MotifParser::new(def_cfg);
    let parser = std::sync::Arc::new(parser);

    // 1.5 Checkpoint resume: if a checkpoint exists, assemble and return (skip scan).
    if let Some(dir) = &config.checkpoint_dir {
        if let Some(existing) = load_checkpoint(dir, signature_def)? {
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
                    Some(scan_one_fragment(fragment, &parser, &signature_def))
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
            pb.set_length(scanned); // total grows as streaming discovers fragments
            pb.inc(chunk_len);
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
    signature_def: &AnnotatorSignatureDefinition,
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
        let ids: Vec<u8> = list.motifs.iter().map(|m| m.id).collect();
        if !signature_def.has_signature(&ids) {
            continue;
        }
        let (frame, strand) = parse_frame(&pseq.identifier);
        for motif in list.motifs {
            let mut m = motif;
            // Java: motif.setDNA(id, offset, seq.getLength(), frame, forwardStrand)
            // Use actual fragment length (not fixed fragment_length).
            m.set_dna(id.clone(), offset, fragment_len, frame, strand);
            out.push((id.clone(), m));
        }
    }
    out
}

/// Resolve mot/store sources (user path takes precedence over built-in) and load MotifDefinition.
fn load_motif_definition(config: &RunConfig) -> std::io::Result<nlr_config::MotifDefinition> {
    let mot_owned;
    let mot_str: &str = match &config.mot_file {
        Some(p) => {
            mot_owned = std::fs::read_to_string(p)
                .map_err(|e| std::io::Error::other(format!("cannot read mot file {}: {}", p.display(), e)))?;
            &mot_owned
        }
        None => EMBEDDED_MOT,
    };
    let store_owned;
    let store_str: &str = match &config.store_file {
        Some(p) => {
            store_owned = std::fs::read_to_string(p)
                .map_err(|e| std::io::Error::other(format!("cannot read store file {}: {}", p.display(), e)))?;
            &store_owned
        }
        None => EMBEDDED_STORE,
    };
    nlr_config::MotifDefinition::load_from_str(mot_str, store_str)
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

/// Remove adjacent duplicates (equivalent to Java `removeRedundantMotifs`: same id + same dnaStart).
///
/// Faithful to Java semantics: each iteration `index++` unconditionally; after a deletion the
/// same position is NOT retried. (Differs from "decrement i after delete": given three identical
/// in a row, Java keeps 2, the decrement variant would reduce to 1.)
fn dedup_adjacent(motifs: &mut Vec<Motif>) {
    let mut i = 0;
    while i < motifs.len() {
        if i + 1 < motifs.len()
            && motifs[i].id == motifs[i + 1].id
            && motifs[i].dna_start == motifs[i + 1].dna_start
        {
            motifs.remove(i + 1);
        }
        i += 1;
    }
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

/// Load checkpoint: read motif results from the checkpoint directory (equivalent to `-c` import).
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
