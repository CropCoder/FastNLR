//! FastNLR command-line entry point.
//!
//! Provides short options `-i/-x/-y/-o/-g/-b/-m/-a/-f/-c/-t/-n`,
//! plus enhancements `-p/--output-prefix/--tmpdir/--progress/--stats/--plot/--summary/--log-level`.
//!
//! `-x`/`-y` are optional: when omitted, built-in mot.txt/store.txt (embedded at compile time)
//! are used, so the tool runs out of the box for standard motif configs.
//!
//! This file owns: arg parsing & help text, input/output validation,
//! pipeline orchestration, startup config summary and end-of-run summary,
//! multi-format output writing.

use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use nlr_cli::{all_motifs, run_with_progress, RunConfig};

/// Valid log-level values (for --log-level validation and hints).
const VALID_LOG_LEVELS: [&str; 5] = ["trace", "debug", "info", "warn", "error"];

/// CLI banner and project metadata shown at the top of `fastnlr --help`.
const HELP_BANNER: &str = r#"   ____           __   _  __ __    ___ 
  / __/___ _ ___ / /_ / |/ // /   / _ \
 / _/ / _ `/(_-</ __//    // /__ / , _/
/_/   \_,_//___/\__//_/|_//____//_/|_| 

FastNLR v1.2.0

Author:     Jiwen Zhao (https://github.com/CropCoder)
Repository: https://github.com/CropCoder/FastNLR
Releases:   https://github.com/CropCoder/FastNLR/releases
Issues:     https://github.com/CropCoder/FastNLR/issues
License:    GPL-3.0-only"#;

/// Print the compact no-argument help block.
fn print_no_args_help() {
    print!("{}", HELP_BANNER);
    print!(
        "\n\nFastNLR: high-speed, accurate NLR immune-receptor locus annotation for plant genomes.\n\n\
         Key parameters:\n  \
         -i INPUT.genomic.fasta   Input genome FASTA (gzip allowed)\n  \
         -o OUTPUT.DIR            Output directory\n  \
         -p PREFIX                Output prefix (default: out)\n  \
         -t THREADS               Thread count\n  \
         --motif-accept-p         Final motif p-value threshold\n  \
         --motif-prelim-p         Motif prefilter p-value threshold\n  \
         --relaxed-seed           Relaxed seeding for divergent NLRs\n\n\
         Usage: fastnlr [OPTIONS] -i <INPUT.genomic.fasta> -o <OUTPUT.DIR>\n\n\
         Use -h for a brief summary or --help for full documentation.\n"
    );
}

#[derive(Parser, Debug)]
#[command(
    name = "fastnlr",
    version,
    before_help = HELP_BANNER,
    long_version = concat!(
        "1.2.0\n",
        "Author:  Jiwen Zhao (https://github.com/CropCoder)\n",
        "Repo:    https://github.com/CropCoder/FastNLR\n",
        "Releases: https://github.com/CropCoder/FastNLR/releases\n",
        "Issues:   https://github.com/CropCoder/FastNLR/issues\n",
        "License:  GPL-3.0-only"
    ),
    about = "FastNLR: high-speed, accurate NLR immune-receptor locus annotation for plant genomes.\n\n\
              Key parameters:\n  \
              -i INPUT.genomic.fasta   Input genome FASTA (gzip allowed)\n  \
              -o OUTPUT.DIR            Output directory\n  \
              -p PREFIX                Output prefix (default: out)\n  \
              -t THREADS               Thread count\n  \
              --motif-accept-p         Final motif p-value threshold\n  \
              --motif-prelim-p         Motif prefilter p-value threshold\n  \
              --relaxed-seed           Relaxed seeding for divergent NLRs",
    long_about = "FastNLR scans six-frame translations of an input FASTA for amino-acid motifs,\n\
                  assembles complete and partial NLR loci, and writes multi-format outputs.\n\
                  By default, the full pipeline runs and writes the complete output set.",
    after_long_help = "Examples:\n  \
        # Full default run (writes the complete output set)\n  \
        fastnlr -i genome.fasta -o results\n\n  \
        # Multithreading with explicit prefix\n  \
        fastnlr -i genome.fasta -o results -p run1 -t 8\n\n  \
        # Checkpoint resume\n  \
        fastnlr -i genome.fasta -o results --checkpoint ckpt/\n\n  \
        # Explicit mot/store override\n  \
        fastnlr -i genome.fasta -x custom_mot.txt -y custom_store.txt -o results"
)]
struct Cli {
    // ===== Input (required) =====
    /// Input genome FASTA (may be gzip-compressed)
    #[arg(
        short = 'i',
        value_name = "INPUT.genomic.fasta",
        help_heading = "Input (required)"
    )]
    input: PathBuf,

    // ===== Motif config (optional, built-in by default) =====
    /// mot.txt PWM config file (default: built-in embedded mot.txt)
    #[arg(short = 'x', help_heading = "Motif config (optional, built-in by default)")]
    mot: Option<PathBuf>,

    /// store.txt CDF config file (default: built-in embedded store.txt)
    #[arg(short = 'y', help_heading = "Motif config (optional, built-in by default)")]
    store: Option<PathBuf>,

    // ===== Output directory and prefix =====
    /// Output directory for all generated files
    #[arg(
        short = 'o',
        required = true,
        value_name = "OUTPUT.DIR",
        help_heading = "Output (required)"
    )]
    output_dir: PathBuf,

    /// Output prefix (default: out); the full output set is written by default
    #[arg(
        short = 'p',
        long = "output-prefix",
        default_value = "out",
        value_name = "PREFIX",
        help_heading = "Output (required)"
    )]
    output_prefix: PathBuf,

    /// Output NLR loci (GFF3)
    #[arg(short = 'g', help_heading = "Optional output files")]
    gff: Option<PathBuf>,

    /// Output NLR loci (BED)
    #[arg(short = 'b', help_heading = "Optional output files")]
    bed: Option<PathBuf>,

    // ===== Sequence / export output =====
    /// Output motif intervals (BED)
    #[arg(short = 'm', help_heading = "Optional sequence/export output")]
    motif_bed: Option<PathBuf>,

    /// Output NB-ARC multiple-alignment fasta
    #[arg(short = 'a', help_heading = "Optional sequence/export output")]
    alignment: Option<PathBuf>,

    /// Output loci sequence fasta; requires 3 args in order: genome.fasta out.fasta flanking(bp)
    #[arg(short = 'f', num_args = 3, help_heading = "Optional sequence/export output")]
    loci: Option<Vec<String>>,

    /// Export precomputed motif results (TSV; usable as checkpoint import)
    #[arg(short = 'c', help_heading = "Optional sequence/export output")]
    export: Option<PathBuf>,

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

    // ===== NLR recall enhancements =====
    /// Final motif-accept p-value threshold (default 1e-5).
    /// Increase (e.g. 1e-3) to improve recall for highly divergent NLRs
    /// (helper NLR/RNL such as ADR1 or NRG1), at the cost of more weak hits.
    #[arg(long, default_value_t = 1e-5, help_heading = "NLR recall enhancements")]
    motif_accept_p: f64,

    /// Motif prefilter p-value threshold (default 1e-4); use together with --motif-accept-p.
    #[arg(long, default_value_t = 1e-4, help_heading = "NLR recall enhancements")]
    motif_prelim_p: f64,

    /// Relaxed seeding: treat a consecutive motif run with at least N NB-ARC motifs
    /// as a seed (disabled by default). Recovers NLR/RNL loci whose motif combinations
    /// are not among the 11 built-in seed combinations.
    #[arg(long, help_heading = "NLR recall enhancements")]
    relaxed_seed: Option<usize>,

    /// Declare domain category for motif IDs outside the built-in tables
    /// (for external mot.txt with extra motifs).
    /// Usage: --motif-category 21=CC,22=CC,25=TIR (repeatable or comma-separated)
    #[arg(long = "motif-category", value_name = "ID=CAT", num_args = 1..)]
    motif_categories: Vec<String>,

    /// Extra seed combinations declared by an external library, e.g. "21,4"
    /// (repeatable; separate combinations with semicolons).
    #[arg(long = "seed-combination", value_name = "IDS", num_args = 1..)]
    extra_seeds: Vec<String>,

    /// Extra signatures (prefilter) declared by an external library, e.g. "21,4".
    #[arg(long = "signature", value_name = "IDS", num_args = 1..)]
    extra_signatures: Vec<String>,
}

fn main() {
    if std::env::args().len() == 1 {
        print_no_args_help();
        return;
    }

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
    config.motif_accept_p = cli.motif_accept_p;
    config.motif_prelim_p = cli.motif_prelim_p;
    config.assemble.relaxed_seed_min_nbarc = cli.relaxed_seed;
    config.motif_categories = cli.motif_categories.clone();
    config.extra_seeds = cli.extra_seeds.clone();
    config.extra_signatures = cli.extra_signatures.clone();

    // Startup config summary: print key params and expected output file list.
    print_run_config(&cli, &config);

    // Progress bar with a rough fragment-total estimate.
    let estimated_fragments =
        estimate_total_fragments(&cli.input, config.fragment_length, config.overlap);
    let pb = make_progress_bar(&cli.progress, estimated_fragments);

    match run_with_progress(&config, pb.as_ref()) {
        Ok(result) => {
            if let Some(pb) = &pb {
                let total = pb.length().unwrap_or(0);
                pb.set_position(total);
                pb.finish();
                std::thread::sleep(std::time::Duration::from_secs(2));
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
            let peak_rss_mb = read_peak_rss_mb();
            print_run_summary(&cli, &config, &result, &written, start.elapsed(), peak_rss_mb);
            tracing::info!("run complete");
        }
        Err(e) => {
            if let Some(pb) = &pb {
                pb.finish_and_clear();
                eprintln!();
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

/// Estimate the total number of chopped fragments from the input file size.
///
/// This is intentionally rough: FASTA headers/newlines and gzip compression
/// mean the estimate is not exact. The final scan position is normalized to
/// 100% when the run completes.
fn estimate_total_fragments(input: &PathBuf, fragment_length: usize, overlap: usize) -> u64 {
    let step = fragment_length.saturating_sub(overlap).max(1) as u64;
    let file_size = std::fs::metadata(input).map(|m| m.len()).unwrap_or(0);
    let estimated_bases = if input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        file_size.saturating_mul(4)
    } else {
        file_size
    };
    estimated_bases.saturating_div(step).max(1)
}

/// Create a progress bar (auto: bar on TTY, otherwise off).
fn make_progress_bar(mode: &str, total: u64) -> Option<indicatif::ProgressBar> {
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
    let pb = indicatif::ProgressBar::new(total);
    // Show: completed/total, percent, elapsed, and remaining time.
    pb.set_style(
        indicatif::ProgressStyle::with_template(
            "[{bar:40}] {percent}% | elapsed {elapsed} | remaining {eta}",
        )
        .unwrap()
        .progress_chars("=> "),
    );
    // 1-second steady tick: refreshes the bar once per second instead of 10×/s,
    // avoiding redraw overhead on fast scans.
    pb.enable_steady_tick(std::time::Duration::from_secs(1));
    Some(pb)
}

/// Pre-validate input and output paths, motif files, and log-level validity.
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

    // The output directory is created before the scan so output writing
    // does not fail after expensive computation.
    if let Err(e) = std::fs::create_dir_all(&cli.output_dir) {
        return Err(format!(
            "cannot create output directory {}: {}",
            cli.output_dir.display(),
            e
        ));
    }

    Ok(())
}

/// Print a CLI error (uniform format, with usage hint).
fn print_cli_error(msg: &str) {
    eprintln!("error: {}", msg);
    eprintln!("run `fastnlr --help` for full usage and examples.");
}

/// Resolve an output path from an optional explicit file or the
/// output-directory + prefix convention.
fn resolve_output_path(
    explicit: &Option<PathBuf>,
    output_dir: &PathBuf,
    prefix: &PathBuf,
    suffix: &str,
) -> PathBuf {
    explicit.clone().unwrap_or_else(|| {
        let base = output_dir.join(prefix);
        PathBuf::from(format!("{}{}", base.display(), suffix))
    })
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
    let resolve = |explicit: &Option<PathBuf>, suffix: &str| {
        resolve_output_path(explicit, &cli.output_dir, &cli.output_prefix, suffix)
    };
    targets.push(resolve(&None::<PathBuf>, ".nlr.txt").display().to_string());
    targets.push(resolve(&cli.gff, ".nlr.gff").display().to_string());
    targets.push(resolve(&cli.bed, ".nlr.bed").display().to_string());
    targets.push(resolve(&cli.motif_bed, ".motifs.bed").display().to_string());
    targets.push(resolve(&cli.alignment, ".nbarc.fasta").display().to_string());
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

/// Read peak resident set size from `/proc/self/status` (Linux).
fn read_peak_rss_mb() -> f64 {
    let status = match std::fs::read_to_string("/proc/self/status") {
        Ok(s) => s,
        Err(_) => return 0.0,
    };
    for line in status.lines() {
        if let Some(value) = line.strip_prefix("VmHWM:") {
            if let Ok(kb) = value.trim().split_whitespace().next().unwrap_or("0").parse::<f64>() {
                return kb / 1024.0;
            }
        }
    }
    0.0
}

/// End-of-run summary: input, resources, NLR counts, type counts, and software metadata.
fn print_run_summary(
    cli: &Cli,
    config: &RunConfig,
    result: &nlr_cli::RunResult,
    written: &[PathBuf],
    elapsed: std::time::Duration,
    peak_rss_mb: f64,
) {
    let def = &result.def;
    let complete = result
        .nlrs
        .iter()
        .filter(|l| l.is_complete_nlr(def))
        .count();
    let mut type_counts: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
    for list in &result.nlrs {
        *type_counts.entry(list.domain_string(def)).or_default() += 1;
    }
    let mut type_counts: Vec<(String, usize)> = type_counts.into_iter().collect();
    type_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    println!("# FastNLR run summary");
    println!(
        "input file\t{}",
        cli.input
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
    );
    println!("threads\t{}", config.threads);
    println!("peak memory\t{:.2} MB", peak_rss_mb);
    println!("elapsed\t{:.2}s", elapsed.as_secs_f64());
    println!("NLR loci\t{}", result.nlrs.len());
    println!("complete NLR\t{}", complete);
    println!("NLR types (count desc):");
    for (class, count) in type_counts {
        println!("{}\t{}", class, count);
    }
    if written.is_empty() {
        println!("output files\t(none)");
    } else {
        println!("output files:");
        for p in written {
            println!("  {}", p.display());
        }
    }
    println!(
        "FastNLR: high-speed, accurate NLR immune-receptor locus annotation for plant genomes | https://github.com/CropCoder/FastNLR"
    );
}

/// Write all output files, returning the list of successfully written paths (for the summary).
/// On the first write failure, return the error and stop.
fn write_outputs(cli: &Cli, result: &nlr_cli::RunResult) -> std::io::Result<Vec<PathBuf>> {
    let def = &result.def;
    let nlrs = &result.nlrs;
    let mut written: Vec<PathBuf> = Vec::new();

    let resolve = |explicit: &Option<PathBuf>, suffix: &str| {
        resolve_output_path(explicit, &cli.output_dir, &cli.output_prefix, suffix)
    };

    // Loci report txt.
    let p = resolve(&None::<PathBuf>, ".nlr.txt");
    let mut f = std::fs::File::create(&p)
        .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
    nlr_output::write_report_txt(&mut f, nlrs, def)
        .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
    written.push(p);
    // Loci GFF.
    let p = resolve(&cli.gff, ".nlr.gff");
    let mut f = std::fs::File::create(&p)
        .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
    let date = now_date_string();
    nlr_output::write_nlr_gff(&mut f, nlrs, def, &date, false)
        .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
    written.push(p);
    // Loci BED.
    let p = resolve(&cli.bed, ".nlr.bed");
    let mut f = std::fs::File::create(&p)
        .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
    nlr_output::write_nlr_bed(&mut f, nlrs, def)
        .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
    written.push(p);
    // Motif BED.
    let p = resolve(&cli.motif_bed, ".motifs.bed");
    let mut f = std::fs::File::create(&p)
        .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
    let motifs = all_motifs(result);
    nlr_output::write_motif_bed(&mut f, &motifs, def, false)
        .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
    written.push(p);
    // NB-ARC alignment fasta.
    let p = resolve(&cli.alignment, ".nbarc.fasta");
    let mut f = std::fs::File::create(&p)
        .map_err(|e| std::io::Error::other(format!("cannot create {}: {}", p.display(), e)))?;
    nlr_output::write_nbarc_alignment_fasta(&mut f, nlrs, def, true)
        .map_err(|e| std::io::Error::other(format!("write {} failed: {}", p.display(), e)))?;
    written.push(p);
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
