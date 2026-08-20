//! FastNLR command-line entry point.
//!
//! Compatible with the original Java flags `-i/-x/-y/-o/-g/-b/-m/-a/-f/-c/-t/-n`,
//! plus enhancements `--output-prefix/--tmpdir/--progress/--stats/--plot/--summary/--log-level`.
//!
//! `-x`/`-y` are optional: when omitted, built-in mot.txt/store.txt (embedded at compile time)
//! are used, so the tool runs out of the box for standard motif configs.
//!
//! This file owns: arg parsing & help text, input validation, pipeline orchestration,
//! startup config summary and end-of-run summary, multi-format output writing.

use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use nlr_cli::{all_motifs, run_with_progress, RunConfig};

/// Valid log-level values (for --log-level validation and hints).
const VALID_LOG_LEVELS: [&str; 5] = ["trace", "debug", "info", "warn", "error"];

#[derive(Parser, Debug)]
#[command(
    name = "fastnlr",
    version,
    long_version = concat!(
        "1.0.0\n",
        "Author:  Jiwen Zhao (https://github.com/CropCoder)\n",
        "Repo:    https://github.com/CropCoder/FastNLR\n",
        "Releases: https://github.com/CropCoder/FastNLR/releases\n",
        "Issues:   https://github.com/CropCoder/FastNLR/issues\n",
        "License:  GPL-3.0-only\n",
        "Original: NLR-Annotator, Steuernagel et al., Plant Physiology, 2020 (PMID: 32184345, GPL-3.0)"
    ),
    about = "FastNLR: scan six reading frames of genomic sequences for amino-acid motifs and annotate NLR immune-receptor loci",
    long_about = "FastNLR: a high-performance Rust rewrite of NLR-Annotator.\n\
                  It scans six-frame translations of an input FASTA for motifs, assembles NLR loci,\n\
                  and emits txt/GFF/BED/alignment-fasta and other formats.\n\
                  Compatible with the original -i/-x/-y/-o/-g/-b/-m/-a/-f/-c/-t/-n flag semantics.\n\n\
                  Repository: https://github.com/CropCoder/FastNLR\n\
                  Author:     Jiwen Zhao (https://github.com/CropCoder)\n\
                  Citation:   Steuernagel et al., Plant Physiology, 2020 (PMID: 32184345)",
    after_long_help = "Examples:\n  \
        # Basic loci annotation\n  \
        fastnlr -i genome.fasta -o out.txt -g out.gff -b out.bed\n\n  \
        # Prefix-derived subfiles + multithreading + plots (uses built-in mot/store)\n  \
        fastnlr -i genome.fasta -p result -t 8 --plot plots/ --summary\n\n  \
        # Checkpoint resume\n  \
        fastnlr -i genome.fasta -p result --checkpoint ckpt/\n\n  \
        # Explicit mot/store override\n  \
        fastnlr -i genome.fasta -x custom_mot.txt -y custom_store.txt -p result\n\n\
        \n\
        Project:   https://github.com/CropCoder/FastNLR\n\
        Releases:  https://github.com/CropCoder/FastNLR/releases\n\
        Issues:    https://github.com/CropCoder/FastNLR/issues\n\
        Author:    Jiwen Zhao (https://github.com/CropCoder)\n\
        License:   GPL-3.0-only  |  Original: NLR-Annotator (Steuernagel et al., 2020, PMID: 32184345)"
)]
struct Cli {
    // ===== Input (required) =====
    /// Input genome FASTA (may be gzip-compressed)
    #[arg(short = 'i', help_heading = "Input (required)")]
    input: PathBuf,

    // ===== Motif config (optional, built-in by default) =====
    /// mot.txt PWM config file (default: built-in embedded mot.txt)
    #[arg(short = 'x', help_heading = "Motif config (optional, built-in by default)")]
    mot: Option<PathBuf>,

    /// store.txt CDF config file (default: built-in embedded store.txt)
    #[arg(short = 'y', help_heading = "Motif config (optional, built-in by default)")]
    store: Option<PathBuf>,

    // ===== Loci output =====
    /// Output NLR loci report (txt)
    #[arg(short = 'o', help_heading = "Loci output")]
    output: Option<PathBuf>,

    /// Output NLR loci (GFF3)
    #[arg(short = 'g', help_heading = "Loci output")]
    gff: Option<PathBuf>,

    /// Output NLR loci (BED)
    #[arg(short = 'b', help_heading = "Loci output")]
    bed: Option<PathBuf>,

    // ===== Sequence / export output =====
    /// Output motif intervals (BED)
    #[arg(short = 'm', help_heading = "Sequence & export output")]
    motif_bed: Option<PathBuf>,

    /// Output NB-ARC multiple-alignment fasta
    #[arg(short = 'a', help_heading = "Sequence & export output")]
    alignment: Option<PathBuf>,

    /// Output loci sequence fasta; requires 3 args in order: genome.fasta out.fasta flanking(bp)
    #[arg(short = 'f', num_args = 3, help_heading = "Sequence & export output")]
    loci: Option<Vec<String>>,

    /// Export precomputed motif results (TSV; usable as checkpoint import)
    #[arg(short = 'c', help_heading = "Sequence & export output")]
    export: Option<PathBuf>,

    // ===== Output control =====
    /// Output prefix; auto-derives .nlr.txt/.nlr.gff/.nlr.bed/.motifs.bed/.nbarc.fasta subfiles
    #[arg(long, short = 'p', help_heading = "Output control")]
    output_prefix: Option<PathBuf>,

    // ===== Performance tuning =====
    /// Thread count (default: auto-detect core count; -t N overrides)
    #[arg(short = 't', help_heading = "Performance tuning")]
    threads: Option<usize>,

    /// Per-thread batch size (fragments; default 1000; affects memory and load granularity)
    #[arg(short = 'n', default_value_t = 1000, help_heading = "Performance tuning")]
    seqs_per_thread: usize,

    // ===== Run enhancements =====
    /// Checkpoint directory (saves motif results after scan; rerun skips scan)
    #[arg(long, help_heading = "Run enhancements")]
    checkpoint: Option<PathBuf>,

    /// Temp file directory (default: system temp dir)
    #[arg(long, help_heading = "Run enhancements")]
    tmpdir: Option<PathBuf>,

    // ===== Observability =====
    /// Progress bar mode: auto (TTY-enabled) / bar / simple / off
    #[arg(long, default_value = "auto", help_heading = "Observability")]
    progress: String,

    /// Run statistics report (TSV: global/per-chromosome/per-motif)
    #[arg(long, help_heading = "Observability")]
    stats: Option<PathBuf>,

    /// Statistics plot output directory (PNG: motif counts, NLRs per chromosome)
    #[arg(long, help_heading = "Observability")]
    plot: Option<PathBuf>,

    /// Per-chromosome human-readable summary (stdout)
    #[arg(long, help_heading = "Observability")]
    summary: bool,

    /// Log level: trace/debug/info/warn/error
    #[arg(long, default_value = "info", help_heading = "Observability")]
    log_level: String,
}

fn main() {
    let start = Instant::now();
    let cli = Cli::parse();
    init_logging(&cli.log_level);

    // Pre-validate inputs before running (fail gracefully, no panic).
    if let Err(msg) = validate_inputs(&cli) {
        print_cli_error(&msg);
        std::process::exit(1);
    }

    // Configure temp dir (--tmpdir > TMPDIR/TEMP > system temp).
    let _tmpdir = setup_tmpdir(&cli);

    let mut config = RunConfig::new(cli.input.clone(), cli.mot.clone(), cli.store.clone());
    if let Some(t) = cli.threads {
        config.threads = t;
    }
    config.seqs_per_thread = cli.seqs_per_thread;
    config.checkpoint_dir = cli.checkpoint.clone();

    // Startup config summary: print key params and expected output file list.
    print_run_config(&cli, &config);

    // Progress bar.
    let pb = make_progress_bar(&cli.progress);

    match run_with_progress(&config, pb.as_ref()) {
        Ok(result) => {
            if let Some(pb) = &pb {
                pb.finish_with_message("scan complete");
            }
            tracing::info!(
                "scan complete: {} sequences, {} NLR loci",
                result.motifs_by_seq.len(),
                result.nlrs.len()
            );

            // Multi-format output (aggregate errors instead of panicking).
            let written = match write_outputs(&cli, &result) {
                Ok(w) => w,
                Err(e) => {
                    print_cli_error(&format!("output write failed: {}", e));
                    std::process::exit(3);
                }
            };

            if cli.summary {
                write_summary(&result);
            }
            if let Some(p) = &cli.stats {
                if let Err(e) = write_stats_file(p, &result) {
                    tracing::warn!("stats report write failed: {}", e);
                }
            }
            if let Some(dir) = &cli.plot {
                if let Err(e) = write_plots(dir, &result) {
                    tracing::warn!("plot generation failed: {}", e);
                }
            }

            // End-of-run summary.
            print_run_summary(&result, &written, start.elapsed());
            tracing::info!("run complete");
        }
        Err(e) => {
            if let Some(pb) = &pb {
                pb.finish_with_message("run failed");
            }
            print_cli_error(&format!("run failed: {}", e));
            std::process::exit(2);
        }
    }
}

/// Initialize logging (invalid level falls back to info).
fn init_logging(level: &str) {
    use tracing_subscriber::EnvFilter;
    let filter = EnvFilter::try_new(level).unwrap_or_else(|_| EnvFilter::new("info"));
    tracing_subscriber::fmt().with_env_filter(filter).init();
}

/// Current UTC time as "YYYY-MM-DD HH:MM:SS" (for the GFF `##date` header).
fn now_date_string() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86400) as i64;
    let rem = secs % 86400;
    let h = rem / 3600;
    let m = (rem % 3600) / 60;
    let s = rem % 60;
    let (y, mo, d) = civil_from_days(days);
    format!("{:04}-{:02}-{:02} {:02}:{:02}:{:02}", y, mo, d, h, m, s)
}

/// Convert days since the Unix epoch (1970-01-01) to a proleptic Gregorian (year, month, day).
/// Howard Hinnant's algorithm (UTC).
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = z - era * 146097; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; // [0, 399]
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = doy - (153 * mp + 2) / 5 + 1; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 }; // [1, 12]
    (if m <= 2 { y + 1 } else { y }, m as u32, d as u32)
}

/// Configure the temp dir and return its path (mainly in-memory here; reserved for checkpoint).
fn setup_tmpdir(cli: &Cli) -> Option<PathBuf> {
    if let Some(dir) = &cli.tmpdir {
        std::fs::create_dir_all(dir).ok()?;
        tracing::info!("temp dir: {}", dir.display());
        return Some(dir.clone());
    }
    if let Ok(dir) = std::env::var("TMPDIR") {
        let p = PathBuf::from(dir);
        std::fs::create_dir_all(&p).ok();
        return Some(p);
    }
    None
}

/// Create a progress bar (auto: bar on TTY, otherwise off).
fn make_progress_bar(mode: &str) -> Option<indicatif::ProgressBar> {
    let enabled = match mode {
        "off" => false,
        "bar" => true,
        "simple" => false, // simple mode degrades to logging
        _ => {
            // auto: enable only on TTY.
            use std::io::IsTerminal;
            std::io::stdout().is_terminal()
        }
    };
    if !enabled {
        return None;
    }
    let pb = indicatif::ProgressBar::new(0);
    // Show: processed/total fragments + percent + throughput (fragments/s) + elapsed + ETA.
    pb.set_style(
        indicatif::ProgressStyle::with_template(
            "{spinner} fragments: {pos}/{len} [{bar:40}] {percent}% | {per_sec} | elapsed {elapsed} eta {eta}",
        )
        .unwrap()
        .progress_chars("=> "),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(100));
    Some(pb)
}

/// Pre-validate: input file existence, log-level validity, no-output-target detection.
fn validate_inputs(cli: &Cli) -> Result<(), String> {
    // Input genome must exist and be a file.
    if !cli.input.is_file() {
        return Err(format!(
            "input file does not exist or is not a file: -i specifies {}",
            cli.input.display()
        ));
    }

    // Optional mot/store: validate only when provided.
    for (flag, path) in [("-x", &cli.mot), ("-y", &cli.store)] {
        if let Some(p) = path {
            if !p.is_file() {
                return Err(format!(
                    "motif config file does not exist or is not a file: {} specifies {}",
                    flag,
                    p.display()
                ));
            }
        }
    }

    // log-level validity.
    if !VALID_LOG_LEVELS.contains(&cli.log_level.as_str()) {
        return Err(format!(
            "invalid --log-level value: {} (valid: {})",
            cli.log_level,
            VALID_LOG_LEVELS.join("/")
        ));
    }

    // No output target -> warn (non-blocking, info only).
    let has_output = cli.output.is_some()
        || cli.gff.is_some()
        || cli.bed.is_some()
        || cli.motif_bed.is_some()
        || cli.alignment.is_some()
        || cli.loci.is_some()
        || cli.export.is_some()
        || cli.output_prefix.is_some()
        || cli.stats.is_some()
        || cli.plot.is_some()
        || cli.summary;
    if !has_output {
        tracing::warn!("no output flags given; will only complete the scan and produce no files");
    }

    Ok(())
}

/// Print a CLI error (uniform format, with usage hint).
fn print_cli_error(msg: &str) {
    eprintln!("error: {}", msg);
    eprintln!("run `fastnlr --help` for full usage and examples.");
}

/// Startup config summary: print key params and expected output file list.
fn print_run_config(cli: &Cli, config: &RunConfig) {
    let mot_src = cli
        .mot
        .as_ref()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|| "built-in".to_string());
    let store_src = cli
        .store
        .as_ref()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|| "built-in".to_string());
    tracing::info!(
        "config: input={} | mot={} | store={} | threads={} | batch={}",
        cli.input.display(),
        mot_src,
        store_src,
        config.threads,
        config.seqs_per_thread
    );

    // Expected output file list.
    let mut targets: Vec<String> = Vec::new();
    let prefix = &cli.output_prefix;
    let resolve = |explicit: &Option<PathBuf>, suffix: &str| -> Option<PathBuf> {
        explicit.clone().or_else(|| {
            prefix
                .as_ref()
                .map(|p| PathBuf::from(format!("{}{}", p.display(), suffix)))
        })
    };
    if let Some(p) = resolve(&cli.output, ".nlr.txt") {
        targets.push(p.display().to_string());
    }
    if let Some(p) = resolve(&cli.gff, ".nlr.gff") {
        targets.push(p.display().to_string());
    }
    if let Some(p) = resolve(&cli.bed, ".nlr.bed") {
        targets.push(p.display().to_string());
    }
    if let Some(p) = resolve(&cli.motif_bed, ".motifs.bed") {
        targets.push(p.display().to_string());
    }
    if let Some(p) = resolve(&cli.alignment, ".nbarc.fasta") {
        targets.push(p.display().to_string());
    }
    if let Some(p) = &cli.export {
        targets.push(p.display().to_string());
    }
    if let Some(p) = &cli.stats {
        targets.push(p.display().to_string());
    }
    if let Some(dir) = &cli.plot {
        targets.push(format!("{}/01-motif-counts.png", dir.display()));
        targets.push(format!("{}/02-chromosome-nlrs.png", dir.display()));
    }
    if !targets.is_empty() {
        tracing::info!("expected outputs: {}", targets.join(", "));
    }
}

/// End-of-run summary: human-readable sequence/NLR/complete counts, elapsed time, output files.
fn print_run_summary(result: &nlr_cli::RunResult, written: &[PathBuf], elapsed: std::time::Duration) {
    let def = &result.def;
    let complete = result
        .nlrs
        .iter()
        .filter(|l| l.is_complete_nlr(def))
        .count();
    println!("# FastNLR done");
    println!("sequences\t{}", result.motifs_by_seq.len());
    println!("NLR loci\t{}", result.nlrs.len());
    println!("complete NLR\t{}", complete);
    println!("elapsed\t{:.2}s", elapsed.as_secs_f64());
    if written.is_empty() {
        println!("output files\t(none)");
    } else {
        println!("output files:");
        for p in written {
            println!("  {}", p.display());
        }
    }
}

/// Write all output files, returning the list of successfully written paths (for the summary).
/// On the first write failure, return the error and stop.
fn write_outputs(cli: &Cli, result: &nlr_cli::RunResult) -> std::io::Result<Vec<PathBuf>> {
    let def = &result.def;
    let nlrs = &result.nlrs;
    let mut written: Vec<PathBuf> = Vec::new();

    let prefix = cli.output_prefix.clone();
    let resolve = |explicit: &Option<PathBuf>, suffix: &str| -> Option<PathBuf> {
        explicit.clone().or_else(|| {
            prefix
                .as_ref()
                .map(|p| PathBuf::from(format!("{}{}", p.display(), suffix)))
        })
    };

    // Loci report txt.
    if let Some(p) = resolve(&cli.output, ".nlr.txt") {
        let mut f = std::fs::File::create(&p)
            .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
        nlr_output::write_report_txt(&mut f, nlrs, def)
            .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
        written.push(p);
    }
    // Loci GFF.
    if let Some(p) = resolve(&cli.gff, ".nlr.gff") {
        let mut f = std::fs::File::create(&p)
            .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
        let date = now_date_string();
        nlr_output::write_nlr_gff(&mut f, nlrs, def, &date, false)
            .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
        written.push(p);
    }
    // Loci BED.
    if let Some(p) = resolve(&cli.bed, ".nlr.bed") {
        let mut f = std::fs::File::create(&p)
            .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
        nlr_output::write_nlr_bed(&mut f, nlrs, def)
            .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
        written.push(p);
    }
    // Motif BED.
    if let Some(p) = resolve(&cli.motif_bed, ".motifs.bed") {
        let mut f = std::fs::File::create(&p)
            .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
        let motifs = all_motifs(result);
        nlr_output::write_motif_bed(&mut f, &motifs, def, false)
            .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
        written.push(p);
    }
    // NB-ARC alignment fasta.
    if let Some(p) = resolve(&cli.alignment, ".nbarc.fasta") {
        let mut f = std::fs::File::create(&p)
            .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
        nlr_output::write_nbarc_alignment_fasta(&mut f, nlrs, def, true)
            .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
        written.push(p);
    }
    // Loci sequence fasta (-f: genome.fasta out.fasta flanking).
    // Extracts from ALL contigs (multi-contig fix), grouping NLRs by contig name.
    if let Some(args) = &cli.loci {
        if args.len() == 3 {
            let genome_path = PathBuf::from(&args[0]);
            let out_path = PathBuf::from(&args[1]);
            let flanking: u64 = match args[2].parse() {
                Ok(v) => v,
                Err(_) => {
                    tracing::warn!("invalid flanking length {}, falling back to 0", args[2]);
                    0
                }
            };
            match nlr_seq::fasta::read_all(&genome_path) {
                Ok(contigs) => {
                    let pairs: Vec<(&str, &str)> = contigs
                        .iter()
                        .map(|c| (c.identifier.as_str(), c.sequence.as_str()))
                        .collect();
                    let mut f = std::fs::File::create(&out_path).map_err(|e| {
                        std::io::Error::other(format!("cannot create {}: {}", out_path.display(), e))
                    })?;
                    nlr_output::write_nlr_loci_all(&mut f, nlrs, &pairs, flanking).map_err(|e| {
                        std::io::Error::other(format!("write {} failed: {}", out_path.display(), e))
                    })?;
                    written.push(out_path);
                }
                Err(e) => {
                    tracing::warn!(
                        "cannot read genome {}, skipping loci fasta output: {}",
                        genome_path.display(),
                        e
                    );
                }
            }
        }
    }
    // Export precomputed motif TSV.
    if let Some(p) = &cli.export {
        let mut f = std::fs::File::create(p)
            .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
        let motifs = all_motifs(result);
        nlr_output::export_motifs(&mut f, &motifs)
            .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
        written.push(p.clone());
    }

    Ok(written)
}

/// Write a statistics report (TSV) via nlr-report.
fn write_stats_file(path: &PathBuf, result: &nlr_cli::RunResult) -> std::io::Result<()> {
    let stats = nlr_report::collect(&result.motifs_by_seq, &result.nlrs, &result.def);
    let mut f = std::fs::File::create(path)
        .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", path.display(), e)))?;
    nlr_report::write_tsv(&mut f, &stats)
        .map_err(|e| std::io::Error::other(format!("stats report write failed: {}", e)))?;
    Ok(())
}

/// Generate statistics plots (motif counts + NLRs per chromosome).
fn write_plots(dir: &PathBuf, result: &nlr_cli::RunResult) -> std::io::Result<()> {
    std::fs::create_dir_all(dir)
        .map_err(|e| std::io::Error::other(format!("cannot create plot dir {}: {}", dir.display(), e)))?;
    let stats = nlr_report::collect(&result.motifs_by_seq, &result.nlrs, &result.def);
    let p1 = dir.join("01-motif-counts.png");
    let p2 = dir.join("02-chromosome-nlrs.png");
    if let Err(e) = nlr_plot::plot_motif_counts(&p1, &stats) {
        tracing::warn!("motif plot generation failed: {}", e);
    }
    if let Err(e) = nlr_plot::plot_chromosome_nlrs(&p2, &stats) {
        tracing::warn!("chromosome plot generation failed: {}", e);
    }
    Ok(())
}

fn write_summary(result: &nlr_cli::RunResult) {
    let def = &result.def;
    println!("#chromosome\tmotifs\tnlrs\tcomplete");
    let mut seq_ids: Vec<&String> = result.motifs_by_seq.keys().collect();
    seq_ids.sort();
    for seq in seq_ids {
        let motif_count = result.motifs_by_seq.get(seq).map(|v| v.len()).unwrap_or(0);
        let nlrs_for_seq: Vec<&nlr_core::motif_list::MotifList> = result
            .nlrs
            .iter()
            .filter(|l| l.sequence_name() == seq.as_str())
            .collect();
        let complete = nlrs_for_seq.iter().filter(|l| l.is_complete_nlr(def)).count();
        println!("{}\t{}\t{}\t{}", seq, motif_count, nlrs_for_seq.len(), complete);
    }
}
