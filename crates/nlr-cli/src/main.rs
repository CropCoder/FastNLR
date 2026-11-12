//! FastNLR command-line entry point.
//!
//! Short options: `-i/-x/-y/-o/-t/-n`. Everything else — `--flank`, `--gff`,
//! `--hmm`, `--checkpoint`, `--progress`, `--log-level`, the recall
//! enhancements and the custom-motif-library declarations — is long-form only;
//! run `fastnlr --help` for the full reference.
//!
//! `-x`/`-y` are optional: when omitted, built-in mot.txt/store.txt (embedded at compile time)
//! are used, so the tool runs out of the box for standard motif configs.
//!
//! This file owns: arg parsing & help text, input/output validation,
//! pipeline orchestration, startup config summary and end-of-run summary,
//! multi-format output writing.

use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Parser;
use nlr_cli::{run_with_progress, AnnotationConfig, RunConfig};

/// Valid log-level values (for --log-level validation and hints).
const VALID_LOG_LEVELS: [&str; 5] = ["trace", "debug", "info", "warn", "error"];

/// 位点序列提取的默认侧翼长度（`--flank` 未给出时）。
const DEFAULT_FLANK: u64 = 2000;

/// CLI banner and project metadata shown at the top of `fastnlr --help`.
/// 版本号取自 Cargo 包版本，避免升级版本时在多处遗漏。
const HELP_BANNER: &str = concat!(
    r#"   ____           __   _  __ __    ___ 
  / __/___ _ ___ / /_ / |/ // /   / _ \
 / _/ / _ `/(_-</ __//    // /__ / , _/
/_/   \_,_//___/\__//_/|_//____//_/|_| 

FastNLR v"#,
    env!("CARGO_PKG_VERSION"),
    r#"

Author:     Jiwen Zhao (https://github.com/CropCoder)
Repository: https://github.com/CropCoder/FastNLR
Releases:   https://github.com/CropCoder/FastNLR/releases
Issues:     https://github.com/CropCoder/FastNLR/issues
License:    GPL-3.0-only"#
);

/// Print the compact no-argument help block.
fn print_no_args_help() {
    print!("{}", HELP_BANNER);
    print!(
        "\n\nFastNLR: high-speed, accurate NLR immune-receptor locus annotation for plant genomes.\n\n\
         Key parameters:\n  \
         -i INPUT.genomic.fasta   Input genome FASTA (gzip allowed)\n  \
         --gff GFF3               Optional GFF3 annotation for gene-protein mode\n  \
         -o OUTPUT.DIR            Output directory\n  \
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
        env!("CARGO_PKG_VERSION"),
        "\n",
        "Author:  Jiwen Zhao (https://github.com/CropCoder)\n",
        "Repo:    https://github.com/CropCoder/FastNLR\n",
        "Releases: https://github.com/CropCoder/FastNLR/releases\n",
        "Issues:   https://github.com/CropCoder/FastNLR/issues\n",
        "License:  GPL-3.0-only"
    ),
    about = "FastNLR: high-speed, accurate NLR immune-receptor locus annotation for plant genomes.\n\n\
              Key parameters:\n  \
              -i INPUT.genomic.fasta   Input genome FASTA (gzip allowed)\n  \
              --gff GFF3               Optional GFF3 annotation for gene-protein mode\n  \
              -o OUTPUT.DIR            Output directory\n  \
              -t THREADS               Thread count\n  \
              --motif-accept-p         Final motif p-value threshold\n  \
              --motif-prelim-p         Motif prefilter p-value threshold\n  \
              --relaxed-seed           Relaxed seeding for divergent NLRs",
    long_about = "FastNLR scans six-frame translations of an input FASTA for amino-acid motifs,\n\
                  assembles complete and partial NLR loci, and writes multi-format outputs.\n\
                  With --gff, FastNLR instead translates the longest CDS transcript per gene\n\
                  and classifies the resulting proteins using ordered HMM domain evidence only.",
    after_long_help = "Examples:\n  \
        # Full default run (writes the complete output set)\n  \
        fastnlr -i genome.fasta -o results\n\n  \
  # GFF3-guided protein annotation (no six-frame translation)\n  \
  fastnlr -i genome.fasta --gff annotation.gff3 -o annotation_results\n\n  \
        # Multithreading\n  \
        fastnlr -i genome.fasta -o results -t 8\n\n  \
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

    /// Optional GFF3 annotation. When supplied, FastNLR translates the longest CDS
    /// transcript per gene instead of performing six-frame genome scanning.
    #[arg(long, value_name = "GFF3", help_heading = "Input (required)")]
    gff: Option<PathBuf>,

    // ===== Motif config (optional, built-in by default) =====
    /// mot.txt PWM config file (default: built-in embedded mot.txt)
    #[arg(
        short = 'x',
        help_heading = "Motif config (optional, built-in by default)"
    )]
    mot: Option<PathBuf>,

    /// store.txt CDF config file (default: built-in embedded store.txt)
    #[arg(
        short = 'y',
        help_heading = "Motif config (optional, built-in by default)"
    )]
    store: Option<PathBuf>,

    // ===== Output directory =====
    /// Output directory for all generated files
    #[arg(
        short = 'o',
        required = true,
        value_name = "OUTPUT.DIR",
        help_heading = "Output (required)"
    )]
    output_dir: PathBuf,

    /// Flanking bases for loci sequence extraction (default 2000; genome mode only)
    #[arg(long, value_name = "BP", help_heading = "Output")]
    flank: Option<u64>,

    // ===== GFF3 annotation mode =====
    /// Extra HMM domain models for GFF3 mode (repeatable, may hold multiple models per file).
    /// A model with the same accession replaces the built-in one; new accessions are added.
    #[arg(long = "hmm", value_name = "FILE", num_args = 1.., help_heading = "GFF3 annotation mode")]
    hmm: Vec<PathBuf>,

    // ===== Performance tuning =====
    /// Thread count (default: auto-detect core count; -t N overrides)
    #[arg(short = 't', help_heading = "Performance tuning")]
    threads: Option<usize>,

    /// Per-thread batch size (fragments; default 1000; affects memory and load granularity)
    #[arg(
        short = 'n',
        default_value_t = 1000,
        help_heading = "Performance tuning"
    )]
    seqs_per_thread: usize,

    // ===== Run enhancements =====
    /// Checkpoint directory (saves motif results after scan; rerun skips scan)
    #[arg(long, help_heading = "Run enhancements")]
    checkpoint: Option<PathBuf>,

    // ===== Observability =====
    /// Progress bar mode: auto (TTY-enabled) / bar / simple / off
    #[arg(long, default_value = "auto", help_heading = "Observability")]
    progress: String,

    /// Log level: trace/debug/info/warn/error
    #[arg(long, default_value = "info", help_heading = "Observability")]
    log_level: String,

    // ===== NLR recall enhancements =====
    /// Final PWM motif-accept p-value threshold (default 1e-5; genome mode only).
    /// Increase (e.g. 1e-3) to improve recall for highly divergent NLRs
    /// (helper NLR/RNL such as ADR1 or NRG1), at the cost of more weak hits.
    #[arg(long, help_heading = "NLR recall enhancements")]
    motif_accept_p: Option<f64>,

    /// PWM motif prefilter p-value threshold (default 1e-4; genome mode only).
    #[arg(long, help_heading = "NLR recall enhancements")]
    motif_prelim_p: Option<f64>,

    /// TIR 专属显著性阈值（默认 1e-7，基因组模式）：剔除 p >= 该值的 TIR 类 motif。
    #[arg(long, value_name = "P", help_heading = "NLR recall enhancements")]
    tir_accept_p: Option<f64>,

    /// Relaxed seeding: treat a consecutive motif run with at least N NB-ARC motifs
    /// as a seed (disabled by default). Recovers NLR/RNL loci whose motif combinations
    /// are not among the built-in exact seed combinations. Mutually exclusive with the
    /// default flexible seeding — pass --no-flexible-seed together with it.
    #[arg(long, help_heading = "NLR recall enhancements")]
    relaxed_seed: Option<usize>,

    /// Flexible seeding threshold: seed when at least N NB-ARC motifs appear in
    /// strictly increasing rank order (a repeated motif does not qualify). Enabled
    /// by default (N=3). Additive to the exact motif-id seed combinations, so loci
    /// found without this flag are never lost.
    #[arg(long, value_name = "N", help_heading = "NLR recall enhancements")]
    flexible_seed: Option<usize>,

    /// Disable flexible seeding and fall back to the exact motif-id seed combinations only.
    #[arg(long, help_heading = "NLR recall enhancements")]
    no_flexible_seed: bool,

    /// Require the P-loop motif (motif_1) for flexible seeding (higher precision).
    #[arg(long, help_heading = "NLR recall enhancements")]
    require_ploop: bool,

    /// Disable high-precision locus filtering and keep all partial loci (legacy behavior).
    #[arg(long, help_heading = "NLR recall enhancements")]
    keep_partial_loci: bool,

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

    if cli.gff.is_some() {
        run_annotation_mode(&cli, start);
        return;
    }

    let mut config = RunConfig::new(cli.input.clone(), cli.mot.clone(), cli.store.clone());
    if let Some(t) = cli.threads {
        config.threads = t;
    }
    config.seqs_per_thread = cli.seqs_per_thread;
    config.checkpoint_dir = cli.checkpoint.clone();
    config.motif_accept_p = cli.motif_accept_p.unwrap_or(1e-5);
    config.motif_prelim_p = cli.motif_prelim_p.unwrap_or(1e-4);
    if let Some(p) = cli.tir_accept_p {
        config.tir_accept_p = p;
    }
    config.assemble.relaxed_seed_min_nbarc = cli.relaxed_seed;
    // 柔性播种默认开启（min_nbarc = 3）；--no-flexible-seed 可关闭，--flexible-seed N 可调阈值。
    config.flexible_seed =
        nlr_cli::resolve_flexible_seed(cli.flexible_seed, cli.no_flexible_seed, cli.require_ploop);
    config.locus_filter = if cli.keep_partial_loci {
        None
    } else {
        Some(nlr_cli::LocusFilterConfig::high_precision_default())
    };
    config.motif_categories = cli.motif_categories.clone();
    config.extra_seeds = cli.extra_seeds.clone();
    config.extra_signatures = cli.extra_signatures.clone();

    // Startup config summary: print key params and expected output file list.
    print_run_config(&cli, &config);

    // 输出前缀：用户输入基因组的文件名（去扩展名）。
    let prefix = genome_prefix(&cli.input);

    // Progress bar with a rough fragment-total estimate.
    let estimated_fragments =
        estimate_total_fragments(&cli.input, config.fragment_length, config.overlap);
    let pb = make_progress_bar(&cli.progress, estimated_fragments);

    match run_with_progress(&config, pb.as_ref()) {
        Ok(mut result) => {
            if let Some(pb) = &pb {
                let total = pb.length().unwrap_or(0);
                pb.set_position(total);
                pb.finish();
            }
            tracing::info!(
                "scan complete: {} sequences, {} NLR loci",
                result.motifs_by_seq.len(),
                result.nlrs.len()
            );

            // 默认输出全部结果（序列 / 表格 / 图 / 摘要）。
            let written =
                match write_outputs(&cli, &prefix, &mut result, config.locus_filter.as_ref()) {
                    Ok(w) => w,
                    Err(e) => {
                        print_cli_error(&format!("output write failed: {}", e));
                        std::process::exit(3);
                    }
                };

            // End-of-run summary.
            let peak_rss_mb = read_peak_rss_mb();
            let summary = build_summary(
                &cli,
                &config,
                &result,
                &written,
                start.elapsed(),
                peak_rss_mb,
            );
            let summary_path = output_path(&cli.output_dir, &prefix, ".summary.txt");
            if let Err(e) = std::fs::write(&summary_path, &summary) {
                tracing::warn!("summary write failed: {}", e);
            }
            print!("{}", summary);
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

/// GFF3 注释模式入口：基因组 + GFF3 → CDS → 蛋白 → HMMER → 写完整结果集。
fn run_annotation_mode(cli: &Cli, start: Instant) {
    let gff = cli.gff.clone().expect("run_annotation_mode requires --gff");
    let mut config = AnnotationConfig::new(cli.input.clone(), gff);
    if let Some(t) = cli.threads {
        config.threads = t;
    }
    config.batch = cli.seqs_per_thread;
    config.hmm_files = cli.hmm.clone();

    let prefix = genome_prefix(&cli.input);
    tracing::info!(
        "GFF3 annotation mode: genome={} | gff={} | built-in HMM library={} ({} models) | extra HMM files={} | threads={} | output={}",
        cli.input.display(),
        config.gff3.display(),
        nlr_domain::MODEL_LIBRARY_VERSION,
        nlr_domain::EMBEDDED_MODELS.len(),
        cli.hmm.len(),
        config.threads,
        cli.output_dir.display()
    );

    let pb = make_progress_bar(&cli.progress, 0);
    let result = match nlr_cli::run_annotation_tracked(&config, |progress| {
        if let Some(pb) = &pb {
            if pb.length() != Some(progress.estimated_total) {
                pb.set_length(progress.estimated_total);
            }
            pb.set_position(progress.scanned_fragments);
        }
        nlr_cli::RunControl::Continue
    }) {
        Ok(r) => r,
        Err(e) => {
            if let Some(pb) = &pb {
                pb.finish_and_clear();
                eprintln!();
            }
            print_cli_error(&format!("GFF3 annotation run failed: {}", e));
            std::process::exit(2);
        }
    };
    if let Some(pb) = &pb {
        let total = pb.length().unwrap_or(0);
        pb.set_position(total);
        pb.finish();
    }

    let meta = nlr_cli::annotation_meta(&config, start.elapsed(), result.model_labels.len());

    let written = match nlr_cli::write_annotation_outputs(&cli.output_dir, &prefix, &result, &meta)
    {
        Ok(w) => w,
        Err(e) => {
            print_cli_error(&format!("output write failed: {}", e));
            std::process::exit(3);
        }
    };

    // 摘要已写入文件，同时打印到 stdout（与 DNA 模式一致）。
    let summary_path = cli.output_dir.join(format!("{prefix}.gff.summary.txt"));
    if let Ok(summary) = std::fs::read_to_string(&summary_path) {
        print!("{summary}");
    }
    let nlr_total = result.records.iter().filter(|r| r.arch.is_nlr).count();
    tracing::info!(
        "GFF3 annotation complete: {} proteins, {} NLR candidates, {} files written",
        result.records.len(),
        nlr_total,
        written.len()
    );
}

/// 时间戳统一格式：年月日时分秒，不保留小数秒。
const DATE_TIME_FORMAT: &str = "%Y-%m-%d %H:%M:%S";

/// 日志时间戳：取运行主机的本地时区（而非 UTC），格式 `YYYY-MM-DD HH:MM:SS`。
struct LocalTimer;

impl tracing_subscriber::fmt::time::FormatTime for LocalTimer {
    fn format_time(&self, w: &mut tracing_subscriber::fmt::format::Writer<'_>) -> std::fmt::Result {
        write!(w, "{}", chrono::Local::now().format(DATE_TIME_FORMAT))
    }
}

/// Initialize logging (invalid level falls back to info).
fn init_logging(level: &str) {
    use std::io::IsTerminal;
    use tracing_subscriber::EnvFilter;
    let filter = EnvFilter::try_new(level).unwrap_or_else(|_| EnvFilter::new("info"));
    // 仅在真正的终端上着色：stdout 被重定向到文件或管道时不写入 ANSI 转义码，
    // 并遵循 NO_COLOR 约定（https://no-color.org）。日志经常被重定向保存，
    // 带着转义码会污染日志文件、破坏下游解析。
    let ansi = std::io::stdout().is_terminal() && std::env::var_os("NO_COLOR").is_none();
    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_timer(LocalTimer)
        .with_ansi(ansi)
        .init();
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
    use std::io::IsTerminal;
    let (enabled, template) = match mode {
        "off" => (false, ""),
        "bar" => (
            true,
            "[{bar:40}] {percent}% | elapsed {elapsed} | remaining {eta}",
        ),
        "simple" => (true, "{pos}/{len} | elapsed {elapsed}"),
        _ => {
            // auto: enable only on TTY.
            (
                std::io::stdout().is_terminal(),
                "[{bar:40}] {percent}% | elapsed {elapsed} | remaining {eta}",
            )
        }
    };
    if !enabled {
        return None;
    }
    let pb = indicatif::ProgressBar::new(total);
    // Show: completed/total, percent, elapsed, and remaining time.
    pb.set_style(
        indicatif::ProgressStyle::with_template(template)
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
    if !has_fasta_extension(&cli.input, &["fa", "fasta", "fna", "fas", "ffn"]) {
        return Err(format!(
            "-i must be a genome FASTA (.fa/.fasta/.fna/.fas/.ffn, optionally .gz): {}",
            cli.input.display()
        ));
    }
    if let Some(gff) = &cli.gff {
        if !gff.is_file() {
            return Err(format!(
                "GFF3 file does not exist or is not a file: --gff specifies {}",
                gff.display()
            ));
        }
        if !has_fasta_extension(gff, &["gff", "gff3"]) {
            return Err(format!(
                "--gff must be a GFF3 file (.gff/.gff3, optionally .gz): {}",
                gff.display()
            ));
        }
    }

    if cli.threads == Some(0) {
        return Err("--threads must be greater than 0".to_string());
    }
    if cli.seqs_per_thread == 0 {
        return Err("-n/--seqs-per-thread must be greater than 0".to_string());
    }
    let motif_accept_p = cli.motif_accept_p.unwrap_or(1e-5);
    let motif_prelim_p = cli.motif_prelim_p.unwrap_or(1e-4);
    for (name, value) in [
        ("--motif-accept-p", motif_accept_p),
        ("--motif-prelim-p", motif_prelim_p),
    ] {
        if !value.is_finite() || value <= 0.0 || value > 1.0 {
            return Err(format!("{name} must be finite and within (0,1]"));
        }
    }
    if motif_prelim_p < motif_accept_p {
        return Err(format!(
            "--motif-prelim-p ({}) must be >= --motif-accept-p ({})",
            motif_prelim_p, motif_accept_p
        ));
    }
    if let Some(p) = cli.tir_accept_p {
        if !p.is_finite() || p <= 0.0 || p > 1.0 {
            return Err("--tir-accept-p must be finite and within (0,1]".to_string());
        }
    }
    if matches!(cli.flexible_seed, Some(0)) || matches!(cli.relaxed_seed, Some(0)) {
        return Err("seed min_nbarc values must be >= 1".to_string());
    }
    if !matches!(cli.progress.as_str(), "auto" | "bar" | "simple" | "off") {
        return Err(format!(
            "invalid --progress value: {} (valid: auto/bar/simple/off)",
            cli.progress
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

    let gff_mode = cli.gff.is_some();

    // DNA 模式播种参数组合校验。柔性播种默认开启，因此 relaxed-seed 必须显式关闭柔性播种。
    if !gff_mode && cli.no_flexible_seed {
        if cli.flexible_seed.is_some() {
            return Err(
                "--flexible-seed 与 --no-flexible-seed 互斥：前者设定阈值，后者关闭柔性播种"
                    .to_string(),
            );
        }
        if cli.require_ploop {
            return Err(
                "--require-ploop 只在柔性播种开启时有意义，不能与 --no-flexible-seed 同时使用"
                    .to_string(),
            );
        }
    } else if !gff_mode && cli.relaxed_seed.is_some() {
        return Err(
            "--relaxed-seed 与默认开启的柔性播种互斥（播种判定柔性优先，该参数会被忽略）；\
             如需宽松播种请改加 --no-flexible-seed"
                .to_string(),
        );
    }

    // GFF3 模式：HMM 文件必须存在，DNA 专属参数必须明确拒绝。
    if gff_mode {
        for path in &cli.hmm {
            if !path.is_file() {
                return Err(format!(
                    "--hmm file does not exist or is not a file: {}",
                    path.display()
                ));
            }
        }
        let mut conflicts: Vec<&str> = Vec::new();
        if cli.mot.is_some() {
            conflicts.push("-x/--mot");
        }
        if cli.store.is_some() {
            conflicts.push("-y/--store");
        }
        if cli.motif_accept_p.is_some() {
            conflicts.push("--motif-accept-p");
        }
        if cli.motif_prelim_p.is_some() {
            conflicts.push("--motif-prelim-p");
        }
        if cli.flank.is_some() {
            conflicts.push("--flank");
        }
        if cli.checkpoint.is_some() {
            conflicts.push("--checkpoint");
        }
        if cli.relaxed_seed.is_some() {
            conflicts.push("--relaxed-seed");
        }
        if cli.flexible_seed.is_some() {
            conflicts.push("--flexible-seed");
        }
        if cli.require_ploop {
            conflicts.push("--require-ploop");
        }
        if cli.no_flexible_seed {
            conflicts.push("--no-flexible-seed");
        }
        if cli.tir_accept_p.is_some() {
            conflicts.push("--tir-accept-p");
        }
        if cli.keep_partial_loci {
            conflicts.push("--keep-partial-loci");
        }
        if !cli.extra_seeds.is_empty() {
            conflicts.push("--seed-combination");
        }
        if !cli.extra_signatures.is_empty() {
            conflicts.push("--signature");
        }
        if !cli.motif_categories.is_empty() {
            conflicts.push("--motif-category");
        }
        if !conflicts.is_empty() {
            return Err(format!(
                "--gff 模式使用纯 HMM 分类，不支持 PWM/DNA 专属参数: {}",
                conflicts.join(", ")
            ));
        }
    } else if !cli.hmm.is_empty() {
        return Err("--hmm 只在 --gff GFF3 注释模式下有效".to_string());
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

fn has_fasta_extension(path: &Path, allowed: &[&str]) -> bool {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let uncompressed = name.strip_suffix(".gz").unwrap_or(&name);
    Path::new(uncompressed)
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| {
            allowed
                .iter()
                .any(|value| extension.eq_ignore_ascii_case(value))
        })
}

/// Print a CLI error (uniform format, with usage hint).
fn print_cli_error(msg: &str) {
    eprintln!("error: {}", msg);
    eprintln!("run `fastnlr --help` for full usage and examples.");
}

/// 从输入基因组文件名推导输出前缀（去掉 .fa/.fasta/.fna/.fas 与 .gz 扩展名）。
fn genome_prefix(input: &Path) -> String {
    let name = input
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("genome");
    let name = name.strip_suffix(".gz").unwrap_or(name);
    let name = name
        .strip_suffix(".fasta")
        .or_else(|| name.strip_suffix(".fna"))
        .or_else(|| name.strip_suffix(".fas"))
        .or_else(|| name.strip_suffix(".fa"))
        .unwrap_or(name);
    name.to_string()
}

/// 在输出目录下按「前缀 + 后缀」拼接输出文件路径。
fn output_path(dir: &Path, prefix: &str, suffix: &str) -> PathBuf {
    dir.join(format!("{}{}", prefix, suffix))
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
    tracing::info!("tir_accept_p: {}", config.tir_accept_p);

    // 播种模式（柔性播种默认开启，便于从日志确认实际生效档位）。
    let seeding = match config.flexible_seed {
        Some(cfg) if cfg.require_ploop => {
            format!("flexible (min_nbarc={}, require_ploop=true)", cfg.min_nbarc)
        }
        Some(cfg) => format!("flexible (min_nbarc={})", cfg.min_nbarc),
        None => match config.assemble.relaxed_seed_min_nbarc {
            Some(n) => format!("relaxed (min_nbarc={})", n),
            None => "exact seed combinations only".to_string(),
        },
    };
    tracing::info!("seeding: {}", seeding);

    // 预期输出文件列表（默认全部输出）。
    let prefix = genome_prefix(&cli.input);
    let suffixes = [
        ".nlr.txt",
        ".nlr.gff",
        ".nlr.bed",
        ".motifs.bed",
        ".nbarc.fasta",
        ".loci.fasta",
        ".stats.tsv",
        ".summary.json",
        ".summary.txt",
    ];
    let targets: Vec<String> = suffixes
        .iter()
        .map(|s| {
            output_path(&cli.output_dir, &prefix, s)
                .display()
                .to_string()
        })
        .collect();
    tracing::info!("expected outputs: {}", targets.join(", "));
}

/// 峰值常驻内存（MB）。
///
/// - Linux：读 `/proc/self/status` 的 `VmHWM`（单位 kB）；
/// - macOS：`getrusage(RUSAGE_SELF).ru_maxrss`（单位字节）；
/// - 其它平台：返回 `None`，摘要里显示 `n/a`。
///
/// 旧实现只支持 Linux，在 macOS/Windows 上恒返回 0，摘要会打印误导性的 `0.00 MB`。
fn read_peak_rss_mb() -> Option<f64> {
    #[cfg(target_os = "linux")]
    {
        let status = std::fs::read_to_string("/proc/self/status").ok()?;
        for line in status.lines() {
            if let Some(value) = line.strip_prefix("VmHWM:") {
                let kb: f64 = value.split_whitespace().next()?.parse().ok()?;
                return Some(kb / 1024.0);
            }
        }
        None
    }

    #[cfg(target_os = "macos")]
    {
        // SAFETY: getrusage 只写入调用方提供的结构体，失败时返回非 0。
        let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
        let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) };
        if rc == 0 {
            // macOS 的 ru_maxrss 单位是字节（Linux 是 kB，故这条分支不共用换算）。
            Some(usage.ru_maxrss as f64 / (1024.0 * 1024.0))
        } else {
            None
        }
    }

    #[cfg(not(any(target_os = "linux", target_os = "macos")))]
    {
        None
    }
}

/// 摘要中的峰值内存文本：取不到时显示 `n/a`，避免出现误导性的 0。
fn format_peak_rss(peak_mb: Option<f64>) -> String {
    match peak_mb {
        Some(mb) => format!("{:.2} MB", mb),
        None => "n/a".to_string(),
    }
}

/// 生成运行摘要文本（写入 .summary.txt 并打印到 stdout）。
fn build_summary(
    cli: &Cli,
    config: &RunConfig,
    result: &nlr_cli::RunResult,
    written: &[PathBuf],
    elapsed: std::time::Duration,
    peak_rss_mb: Option<f64>,
) -> String {
    use std::fmt::Write as _;
    let def = &result.def;
    let mut type_counts: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    let mut family_counts: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    for list in &result.nlrs {
        *type_counts.entry(list.domain_string(def)).or_default() += 1;
        *family_counts
            .entry(list.family(def).to_string())
            .or_default() += 1;
    }
    let mut type_counts: Vec<(String, usize)> = type_counts.into_iter().collect();
    type_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let mut family_counts: Vec<(String, usize)> = family_counts.into_iter().collect();
    family_counts.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let mut s = String::new();
    let _ = writeln!(s, "# FastNLR run summary");
    let _ = writeln!(
        s,
        "# input file\t{}",
        cli.input
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
    );
    let _ = writeln!(s, "# threads\t{}", config.threads);
    let _ = writeln!(s, "# peak memory\t{}", format_peak_rss(peak_rss_mb));
    let _ = writeln!(s, "# elapsed\t{:.2}s", elapsed.as_secs_f64());
    let _ = writeln!(s, "# NLR loci\t{}", result.nlrs.len());
    let _ = writeln!(
        s,
        "# Complete NLR loci\t{}",
        result
            .nlrs
            .iter()
            .filter(|list| list.is_complete_nlr(def))
            .count()
    );
    let _ = writeln!(s, "# NLR types (count desc):");
    for (class, count) in type_counts {
        let _ = writeln!(s, "# {}\t{}", count, class);
    }
    let _ = writeln!(s, "# NLR families (count desc):");
    for (family, count) in family_counts {
        let _ = writeln!(s, "# {}\t{}", count, family);
    }
    let _ = writeln!(s, "# output files:");
    for p in written {
        let _ = writeln!(s, "# {}", p.display());
    }
    let _ = writeln!(
        s,
        "# FastNLR: high-speed, accurate NLR immune-receptor locus annotation for plant genomes | https://github.com/CropCoder/FastNLR"
    );
    s
}

/// 写入全部输出文件（序列 / 表格），返回成功写入的路径列表。
fn write_outputs(
    cli: &Cli,
    prefix: &str,
    result: &mut nlr_cli::RunResult,
    locus_filter: Option<&nlr_cli::LocusFilterConfig>,
) -> std::io::Result<Vec<PathBuf>> {
    nlr_cli::write_genome_outputs(
        &cli.input,
        &cli.output_dir,
        prefix,
        cli.flank.unwrap_or(DEFAULT_FLANK),
        result,
        locus_filter,
    )
}

#[cfg(test)]
mod tests {
    use super::{format_peak_rss, read_peak_rss_mb};

    #[test]
    fn peak_rss_text_is_explicit_when_unavailable() {
        assert_eq!(format_peak_rss(Some(2048.0)), "2048.00 MB");
        assert_eq!(format_peak_rss(None), "n/a");
    }

    #[test]
    #[cfg(target_os = "linux")]
    fn peak_rss_is_available_on_linux() {
        let v = read_peak_rss_mb();
        assert!(
            v.is_some_and(|mb| mb > 0.0),
            "Linux 上应能读到峰值内存: {v:?}"
        );
    }
}
