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

    // 2. Chop all fragments.
    let chopper = nlr_seq::SequenceChopper::from_file(
        &config.input_fasta,
        config.fragment_length,
        config.overlap,
    )?;
    let mut fragments = Vec::new();
    let mut chopper = chopper;
    while let Some(fragment) = chopper.next_sequence() {
        fragments.push(fragment);
    }
    tracing::info!("chopping done: {} fragments", fragments.len());

    // 3. Parallel batch scan (seqs_per_thread fragments per batch).
    let num_threads = config.threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .map_err(|e| std::io::Error::other(e.to_string()))?;

    let scan_batch = |batch: &[nlr_seq::translate::BioSequence]| -> Vec<(String, Motif)> {
        let mut local: HashMap<String, Vec<Motif>> = HashMap::new();
        for fragment in batch {
            let (id, offset) = parse_fragment_id(&fragment.identifier);
            let fragment_len = fragment.len() as u64; // actual fragment length (last fragment may be < fragment_length)
            let protein_seqs = fragment.translate2protein();
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
                    local.entry(id.clone()).or_default().push(m);
                }
            }
        }
        local.into_iter().flat_map(|(k, v)| v.into_iter().map(move |m| (k.clone(), m))).collect()
    };

    let batches: Vec<&[nlr_seq::translate::BioSequence]> =
        fragments.chunks(config.seqs_per_thread.max(1)).collect();

    // Progress bar granularity = fragments; total length = fragment count.
    if let Some(pb) = progress {
        pb.set_length(fragments.len() as u64);
    }

    let results: Vec<Vec<(String, Motif)>> = pool.install(|| {
        use rayon::prelude::*;
        batches
            .par_iter()
            .filter_map(|batch| {
                // Skip remaining batches after interrupt.
                if interrupted.load(std::sync::atomic::Ordering::SeqCst) {
                    return None;
                }
                let r = scan_batch(batch);
                if let Some(pb) = progress {
                    pb.inc(batch.len() as u64); // advance by batch.len() fragments
                }
                Some(r)
            })
            .collect()
    });

    if interrupted.load(std::sync::atomic::Ordering::SeqCst) {
        tracing::warn!("scan interrupted, outputting completed batch results");
    }

    // 4. Merge + dedup.
    let mut motifs_by_seq: HashMap<String, Vec<Motif>> = HashMap::new();
    for batch in results {
        for (seq_id, motif) in batch {
            motifs_by_seq.entry(seq_id).or_default().push(motif);
        }
    }
    for v in motifs_by_seq.values_mut() {
        v.sort();
        dedup_adjacent(v);
    }

    // 5. Assemble NLRs.
    let mut nlrs: Vec<MotifList> = Vec::new();
    let mut seq_ids: Vec<&String> = motifs_by_seq.keys().collect();
    seq_ids.sort();
    for seq_id in seq_ids {
        let motifs = motifs_by_seq.get(seq_id).unwrap().clone();
        let mut assembled = nlr_assemble::assemble(seq_id, motifs, &config.assemble, &signature_def);
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
