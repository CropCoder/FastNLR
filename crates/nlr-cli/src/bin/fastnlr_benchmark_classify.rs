//! Benchmark-only protein FASTA classifier.
//!
//! This helper reuses the same ordered HMM-only rules as FastNLR. It is intentionally
//! separate from the user-facing CLI so genomic mode remains the only supported input.
//!
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use std::io::Write;
use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Parser)]
#[command(
    name = "fastnlr-benchmark-classify",
    about = "Classify protein FASTA records with FastNLR ordered HMM rules"
)]
struct Args {
    /// Protein FASTA input.
    #[arg(short = 'i', long)]
    input: PathBuf,

    /// Output prefix; writes <prefix>.classification.tsv and <prefix>.domains.tsv.
    #[arg(short = 'o', long)]
    output_prefix: PathBuf,

    /// Worker threads.
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Proteins per batch.
    #[arg(short = 'n', long, default_value_t = 1000)]
    batch: usize,

    /// Additional HMM files, repeatable.
    #[arg(long = "hmm", num_args = 1..)]
    hmm: Vec<PathBuf>,
}

struct Record {
    id: String,
    length: usize,
    arch: nlr_domain::NlrArchitecture,
    hits: Vec<nlr_domain::DomainHit>,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args = Args::parse();
    if args.threads == 0 || args.batch == 0 {
        return Err("--threads and --batch must be greater than zero".to_string());
    }
    let mut library = nlr_domain::HmmLibrary::from_embedded();
    for path in &args.hmm {
        library
            .add_file(path)
            .map_err(|error| format!("cannot load {}: {error}", path.display()))?;
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .map_err(|error| error.to_string())?;
    let mut reader = nlr_seq::fasta::FastaReader::from_file(&args.input)
        .map_err(|error| format!("cannot read {}: {error}", args.input.display()))?;
    let mut records = Vec::new();
    loop {
        let mut batch = Vec::with_capacity(args.batch);
        while batch.len() < args.batch {
            match reader
                .read_entry()
                .map_err(|error| format!("FASTA read failed: {error}"))?
            {
                Some(sequence) => batch.push(sequence),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        let mut chunk: Vec<Record> = pool.install(|| {
            use rayon::prelude::*;
            batch
                .par_iter()
                .map_init(
                    || library.protein_scanner(),
                    |scanner, sequence| {
                        let hits = scanner.scan_protein(&sequence.identifier, &sequence.sequence);
                        let arch = nlr_domain::classify(&hits);
                        Record {
                            id: sequence.identifier.clone(),
                            length: sequence.sequence.len(),
                            arch,
                            hits,
                        }
                    },
                )
                .collect()
        });
        records.append(&mut chunk);
    }
    records.sort_by(|left, right| left.id.cmp(&right.id));

    let prefix = args.output_prefix.display().to_string();
    let classification_path = format!("{prefix}.classification.tsv");
    let domains_path = format!("{prefix}.domains.tsv");
    if let Some(parent) = args.output_prefix.parent() {
        std::fs::create_dir_all(parent)
            .map_err(|error| format!("cannot create {}: {error}", parent.display()))?;
    }
    let mut classification = std::fs::File::create(&classification_path)
        .map_err(|error| format!("cannot create {classification_path}: {error}"))?;
    writeln!(
        classification,
        "#protein_id\tlength\tclass\tfamily\tnterm\thas_nbarc\thas_tir\thas_lrr_any\thas_cterm_lrr\tcomplete\tn_domains\tdomains"
    )
    .map_err(|error| error.to_string())?;
    let mut domains = std::fs::File::create(&domains_path)
        .map_err(|error| format!("cannot create {domains_path}: {error}"))?;
    writeln!(
        domains,
        "#protein_id\tsource\tmodel\tmodel_name\tcategory\tstart\tend\tbitscore\tpvalue\tpassed"
    )
    .map_err(|error| error.to_string())?;

    for record in &records {
        let passed = record
            .hits
            .iter()
            .filter(|hit| hit.passed)
            .map(|hit| {
                (
                    hit.start,
                    format!("{}:{}-{}", hit.model_name, hit.start, hit.end),
                )
            })
            .collect::<Vec<_>>();
        let domain_count = passed.len();
        let mut passed = passed;
        passed.sort_by_key(|(start, _)| *start);
        writeln!(
            classification,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.id,
            record.length,
            record.arch.class,
            record.arch.family,
            record.arch.nterm,
            record.arch.has_nbarc,
            record.arch.has_tir,
            record.arch.has_lrr_any,
            record.arch.has_cterm_lrr,
            record.arch.complete,
            domain_count,
            passed
                .iter()
                .map(|(_, label)| label.as_str())
                .collect::<Vec<_>>()
                .join(";")
        )
        .map_err(|error| error.to_string())?;
        for hit in &record.hits {
            writeln!(
                domains,
                "{}\thmm\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}",
                record.id,
                hit.model_acc,
                hit.model_name,
                hit.category.as_str(),
                hit.start,
                hit.end,
                hit.bitscore,
                nlr_core::motif::format_double_java(hit.pvalue),
                hit.passed
            )
            .map_err(|error| error.to_string())?;
        }
    }
    println!(
        "classified {} proteins -> {}, {}",
        records.len(),
        classification_path,
        domains_path
    );
    Ok(())
}
