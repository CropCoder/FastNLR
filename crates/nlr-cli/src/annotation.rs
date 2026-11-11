//! GFF3 驱动的基因蛋白注释流水线。
//!
//! 该流程不执行六框翻译：先流式解析 GFF3 并选择每个基因的最长转录本，再按 contig
//! 流式读取基因组、提取 CDS、翻译蛋白，最后使用纯 HMM 结构域证据进行有序分类。
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as _;
use std::io::Write;
use std::path::{Path, PathBuf};

use nlr_core::motif::format_double_java;
use nlr_domain::{DomainHit, NlrArchitecture};
use nlr_seq::gff3::{read_gff3, AnnotationSet, SkippedAnnotation, Strand, TranscriptModel};
use serde::Serialize;

use crate::{RunControl, RunProgress};

/// GFF3 注释模式配置。
#[derive(Debug, Clone)]
pub struct AnnotationConfig {
    pub genome_fasta: PathBuf,
    pub gff3: PathBuf,
    pub hmm_files: Vec<PathBuf>,
    pub threads: usize,
    /// 每批翻译和扫描的基因数。
    pub batch: usize,
}

impl AnnotationConfig {
    pub fn new(genome_fasta: PathBuf, gff3: PathBuf) -> Self {
        AnnotationConfig {
            genome_fasta,
            gff3,
            hmm_files: Vec::new(),
            threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1),
            batch: 1000,
        }
    }

    pub fn validate(&self) -> std::io::Result<()> {
        let invalid =
            |message: String| std::io::Error::new(std::io::ErrorKind::InvalidInput, message);
        if !self.genome_fasta.is_file() {
            return Err(invalid(format!(
                "genome FASTA does not exist: {}",
                self.genome_fasta.display()
            )));
        }
        if !self.gff3.is_file() {
            return Err(invalid(format!(
                "GFF3 file does not exist: {}",
                self.gff3.display()
            )));
        }
        if self.threads == 0 {
            return Err(invalid("thread count must be greater than 0".into()));
        }
        if self.batch == 0 {
            return Err(invalid("batch size must be greater than 0".into()));
        }
        Ok(())
    }
}

/// 结构域在基因组上的一个片段；坐标均为 1-based、闭区间。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicPiece {
    pub start: u64,
    pub end: u64,
}

/// HMM 结构域证据。
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationDomain {
    pub source: String,
    pub model: String,
    pub model_name: String,
    pub category: String,
    pub protein_start: usize,
    pub protein_end: usize,
    pub cds_start: usize,
    pub cds_end: usize,
    pub genomic_start: u64,
    pub genomic_end: u64,
    pub genomic_pieces: Vec<GenomicPiece>,
    pub bitscore: Option<f32>,
    pub pvalue: Option<f64>,
    pub passed: bool,
}

/// 单个最长转录本及其注释结果。
#[derive(Debug, Clone)]
pub struct AnnotationRecord {
    pub transcript: TranscriptModel,
    pub cds_sequence: String,
    pub protein_sequence: String,
    pub internal_stop_count: usize,
    pub has_internal_stop: bool,
    pub trailing_partial_bases: u8,
    pub arch: NlrArchitecture,
    pub domains: Vec<AnnotationDomain>,
}

impl AnnotationRecord {
    pub fn protein_length(&self) -> usize {
        self.protein_sequence.len()
    }

    pub fn passed_domain_count(&self) -> usize {
        self.domains.iter().filter(|domain| domain.passed).count()
    }

    pub fn domain_layout(&self) -> String {
        let mut domains: Vec<&AnnotationDomain> =
            self.domains.iter().filter(|domain| domain.passed).collect();
        domains.sort_by_key(|domain| (domain.protein_start, domain.protein_end));
        domains
            .iter()
            .map(|domain| {
                format!(
                    "{}:{}-{}",
                    domain.model, domain.protein_start, domain.protein_end
                )
            })
            .collect::<Vec<_>>()
            .join(";")
    }
}

/// GFF3 模式整体结果。
#[derive(Debug, Clone, Default)]
pub struct AnnotationResult {
    pub records: Vec<AnnotationRecord>,
    pub skipped: Vec<SkippedAnnotation>,
    pub genes_with_cds: usize,
    pub model_labels: Vec<String>,
}

/// 输出摘要所需运行元信息。
#[derive(Debug, Clone)]
pub struct AnnotationRunMeta {
    pub genome_name: String,
    pub gff_name: String,
    pub threads: usize,
    pub elapsed: std::time::Duration,
    pub model_library_version: String,
    pub model_count: usize,
    pub hmm_files: Vec<String>,
}

/// 无回调运行的便捷入口。
pub fn run_annotation(config: &AnnotationConfig) -> std::io::Result<AnnotationResult> {
    run_annotation_tracked(config, |_| RunControl::Continue)
}

/// 流式、并行执行 GFF3 注释；每处理完一个基因批次回调一次。
pub fn run_annotation_tracked(
    config: &AnnotationConfig,
    mut on_progress: impl FnMut(RunProgress) -> RunControl,
) -> std::io::Result<AnnotationResult> {
    config.validate()?;

    let AnnotationSet {
        transcripts,
        mut skipped,
        genes_with_cds,
    } = read_gff3(&config.gff3)?;
    let estimated_total = transcripts.len() as u64;
    tracing::info!(
        "GFF3 parsed: {} genes with CDS, {} longest transcripts",
        genes_with_cds,
        transcripts.len()
    );

    let mut by_seqid: HashMap<String, Vec<TranscriptModel>> = HashMap::new();
    for transcript in transcripts {
        by_seqid
            .entry(transcript.seqid.clone())
            .or_default()
            .push(transcript);
    }
    for models in by_seqid.values_mut() {
        models.sort_by(|a, b| (a.start0, &a.gene.gene_id).cmp(&(b.start0, &b.gene.gene_id)));
    }

    let mut library = nlr_domain::HmmLibrary::from_embedded();
    for path in &config.hmm_files {
        let added = library.add_file(path)?;
        tracing::info!("loaded {added} domain model(s) from {}", path.display());
    }
    let model_labels: Vec<String> = library
        .models()
        .iter()
        .map(|model| format!("{}({})", model.acc, model.name))
        .collect();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads.max(1))
        .build()
        .map_err(|error| std::io::Error::other(error.to_string()))?;
    let interrupted = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
    let interrupted_flag = interrupted.clone();
    let _ = ctrlc::set_handler(move || {
        tracing::warn!("interrupt signal received, stopping GFF3 annotation after current batch");
        interrupted_flag.store(true, std::sync::atomic::Ordering::SeqCst);
    });

    let mut reader = nlr_seq::fasta::FastaReader::from_file(&config.genome_fasta)?;
    let mut records = Vec::with_capacity(estimated_total as usize);
    let mut processed = 0u64;
    let mut canceled = false;
    let mut seen_seqids = HashSet::new();
    let batch_size = config.batch.max(1);

    while let Some(contig) = reader.read_entry()? {
        let Some(models) = by_seqid.remove(&contig.identifier) else {
            continue;
        };
        seen_seqids.insert(contig.identifier.clone());
        for chunk in models.chunks(batch_size) {
            let chunk_results: Vec<Result<AnnotationRecord, (String, String, String, String)>> =
                pool.install(|| {
                    use rayon::prelude::*;
                    chunk
                        .par_iter()
                        .map_init(
                            || library.protein_scanner(),
                            |scanner, transcript| {
                                scan_transcript(scanner, transcript, &contig.sequence).map_err(
                                    |reason| {
                                        (
                                            transcript.gene.gene_id.clone(),
                                            transcript.transcript_id.clone(),
                                            transcript.seqid.clone(),
                                            reason,
                                        )
                                    },
                                )
                            },
                        )
                        .collect()
                });

            for item in chunk_results {
                match item {
                    Ok(record) => records.push(record),
                    Err((gene_id, transcript_id, seqid, reason)) => {
                        skipped.push(SkippedAnnotation {
                            gene_id,
                            transcript_id,
                            seqid,
                            reason,
                        });
                    }
                }
            }
            processed += chunk.len() as u64;
            if matches!(
                on_progress(RunProgress {
                    scanned_fragments: processed,
                    estimated_total,
                }),
                RunControl::Break
            ) || interrupted.load(std::sync::atomic::Ordering::SeqCst)
            {
                canceled = true;
                break;
            }
        }
        if canceled {
            break;
        }
    }

    if !canceled {
        for models in by_seqid.into_values() {
            for transcript in models {
                skipped.push(SkippedAnnotation {
                    gene_id: transcript.gene.gene_id,
                    transcript_id: transcript.transcript_id,
                    seqid: transcript.seqid,
                    reason: "contig_not_found_in_genome_fasta".to_string(),
                });
            }
        }
    }

    records.sort_by(|a, b| {
        (
            &a.transcript.seqid,
            a.transcript.start0,
            &a.transcript.gene.gene_id,
        )
            .cmp(&(
                &b.transcript.seqid,
                b.transcript.start0,
                &b.transcript.gene.gene_id,
            ))
    });
    skipped.sort_by(|a, b| {
        (&a.seqid, &a.gene_id, &a.transcript_id).cmp(&(&b.seqid, &b.gene_id, &b.transcript_id))
    });
    tracing::info!(
        "GFF3 annotation scan complete: {} proteins, {} skipped, {} contigs matched",
        records.len(),
        skipped.len(),
        seen_seqids.len()
    );

    Ok(AnnotationResult {
        records,
        skipped,
        genes_with_cds,
        model_labels,
    })
}

fn scan_transcript(
    scanner: &mut nlr_domain::ProteinScanner,
    transcript: &TranscriptModel,
    contig_sequence: &str,
) -> Result<AnnotationRecord, String> {
    let cds_sequence = transcript
        .extract_cds(contig_sequence)
        .map_err(|error| error.to_string())?;
    let translated = translate_cds(&cds_sequence, transcript.translation_offset());
    if translated.protein.is_empty() {
        return Err("translated protein is empty".to_string());
    }

    let hits = scanner.scan_protein(&transcript.transcript_id, &translated.protein);
    let arch = nlr_domain::classify(&hits);
    let domains = build_domains(transcript, &hits)?;

    Ok(AnnotationRecord {
        transcript: transcript.clone(),
        cds_sequence,
        protein_sequence: translated.protein,
        internal_stop_count: translated.internal_stop_count,
        has_internal_stop: translated.internal_stop_count > 0,
        trailing_partial_bases: translated.trailing_partial_bases,
        arch,
        domains,
    })
}

struct TranslatedCds {
    protein: String,
    internal_stop_count: usize,
    trailing_partial_bases: u8,
}

// `as_chunks` is newer than the workspace's Rust 1.85 minimum supported version.
#[allow(clippy::chunks_exact_to_as_chunks)]
fn translate_cds(cds: &str, translation_offset: usize) -> TranslatedCds {
    let translation_offset = translation_offset.min(cds.len());
    let coding_bases = cds.len() - translation_offset;
    let trailing_partial_bases = (coding_bases % 3) as u8;
    let coding_end = cds.len() - trailing_partial_bases as usize;
    let mut protein = String::with_capacity((coding_end - translation_offset) / 3);
    for codon in cds.as_bytes()[translation_offset..coding_end].chunks_exact(3) {
        protein.push(nlr_seq::codon::translate_triplet(codon));
    }
    if protein.ends_with('*') {
        protein.pop();
    }
    let internal_stop_count = protein
        .chars()
        .filter(|character| *character == '*')
        .count();
    let protein = protein
        .chars()
        .map(|character| if character == '*' { 'X' } else { character })
        .collect();
    TranslatedCds {
        protein,
        internal_stop_count,
        trailing_partial_bases,
    }
}

fn build_domains(
    transcript: &TranscriptModel,
    hits: &[DomainHit],
) -> Result<Vec<AnnotationDomain>, String> {
    let mut domains = Vec::with_capacity(hits.len());
    for hit in hits {
        domains.push(annotation_domain(
            transcript,
            "hmm".to_string(),
            hit.model_acc.clone(),
            hit.model_name.clone(),
            hit.category.as_str().to_string(),
            hit.start,
            hit.end,
            Some(hit.bitscore),
            Some(hit.pvalue),
            hit.passed,
        )?);
    }
    domains.sort_by_key(|domain| (domain.protein_start, domain.protein_end));
    Ok(domains)
}

#[allow(clippy::too_many_arguments)]
fn annotation_domain(
    transcript: &TranscriptModel,
    source: String,
    model: String,
    model_name: String,
    category: String,
    protein_start: usize,
    protein_end: usize,
    bitscore: Option<f32>,
    pvalue: Option<f64>,
    passed: bool,
) -> Result<AnnotationDomain, String> {
    let (cds_start, cds_end, genomic_start, genomic_end, genomic_pieces) =
        map_protein_interval(transcript, protein_start, protein_end)
            .ok_or_else(|| format!("{model} coordinates are outside the translated CDS"))?;
    Ok(AnnotationDomain {
        source,
        model,
        model_name,
        category,
        protein_start,
        protein_end,
        cds_start,
        cds_end,
        genomic_start,
        genomic_end,
        genomic_pieces,
        bitscore,
        pvalue,
        passed,
    })
}

fn map_protein_interval(
    transcript: &TranscriptModel,
    protein_start: usize,
    protein_end: usize,
) -> Option<(usize, usize, u64, u64, Vec<GenomicPiece>)> {
    if protein_start == 0 || protein_end < protein_start {
        return None;
    }
    let cds_len = transcript.cds_length_bp as usize;
    let translation_offset = transcript.translation_offset().min(cds_len);
    let interval_start0 = translation_offset
        .checked_add(protein_start.checked_sub(1)?.checked_mul(3)?)?
        .min(cds_len);
    let interval_end0 = translation_offset
        .checked_add(protein_end.checked_mul(3)?)?
        .min(cds_len);
    if interval_start0 >= interval_end0 {
        return None;
    }

    let mut cursor = 0usize;
    let mut pieces = Vec::new();
    for segment in &transcript.cds {
        let segment_len = usize::try_from(segment.len()).map_err(|_| ()).ok()?;
        let segment_start = cursor;
        let segment_end = cursor + segment_len;
        cursor = segment_end;
        let overlap_start = interval_start0.max(segment_start);
        let overlap_end = interval_end0.min(segment_end);
        if overlap_start >= overlap_end {
            continue;
        }
        let local_start = overlap_start - segment_start;
        let local_end = overlap_end - segment_start;
        let (start0, end0) = match transcript.strand {
            Strand::Plus => (
                segment.start0 + local_start as u64,
                segment.start0 + local_end as u64,
            ),
            Strand::Minus => (
                segment.end0.saturating_sub(local_end as u64),
                segment.end0.saturating_sub(local_start as u64),
            ),
        };
        pieces.push(GenomicPiece {
            start: start0 + 1,
            end: end0,
        });
    }
    if pieces.is_empty() {
        return None;
    }
    let genomic_start = pieces.iter().map(|piece| piece.start).min()?;
    let genomic_end = pieces.iter().map(|piece| piece.end).max()?;
    Some((
        interval_start0 + 1,
        interval_end0,
        genomic_start,
        genomic_end,
        pieces,
    ))
}

/// 写出 GFF3 模式的完整结果集。
pub fn write_annotation_outputs(
    output_dir: &Path,
    prefix: &str,
    result: &AnnotationResult,
    meta: &AnnotationRunMeta,
) -> std::io::Result<Vec<PathBuf>> {
    std::fs::create_dir_all(output_dir)?;
    let mut written = Vec::new();

    let cds_path = output_dir.join(format!("{prefix}.gff.cds.fasta"));
    write_annotation_fasta(
        &cds_path,
        result.records.iter().map(|record| {
            (
                record.transcript.transcript_id.clone(),
                format_record_description(record),
                record.cds_sequence.clone(),
            )
        }),
    )?;
    written.push(cds_path);

    let protein_path = output_dir.join(format!("{prefix}.gff.proteins.fasta"));
    write_annotation_fasta(
        &protein_path,
        result.records.iter().map(|record| {
            (
                record.transcript.transcript_id.clone(),
                format_record_description(record),
                record.protein_sequence.clone(),
            )
        }),
    )?;
    written.push(protein_path);

    let domains_path = output_dir.join(format!("{prefix}.gff.domains.tsv"));
    write_domains_tsv(&domains_path, &result.records)?;
    written.push(domains_path);

    let genes_path = output_dir.join(format!("{prefix}.gff.genes.tsv"));
    write_records_tsv(&genes_path, &result.records)?;
    written.push(genes_path);

    let nlr_path = output_dir.join(format!("{prefix}.gff.nlr.tsv"));
    let nlr_records: Vec<&AnnotationRecord> = result
        .records
        .iter()
        .filter(|record| record.arch.is_nlr)
        .collect();
    write_record_refs_tsv(&nlr_path, &nlr_records)?;
    written.push(nlr_path);

    let nlr_gff_path = output_dir.join(format!("{prefix}.gff.nlr.gff3"));
    write_nlr_gff3(&nlr_gff_path, &nlr_records)?;
    written.push(nlr_gff_path);

    let skipped_path = output_dir.join(format!("{prefix}.gff.skipped.tsv"));
    write_skipped_tsv(&skipped_path, &result.skipped)?;
    written.push(skipped_path);

    let summary_json_path = output_dir.join(format!("{prefix}.gff.summary.json"));
    write_summary_json(&summary_json_path, result, meta)?;
    written.push(summary_json_path);

    let summary_path = output_dir.join(format!("{prefix}.gff.summary.txt"));
    write_summary_text(&summary_path, &written, result, meta)?;
    written.push(summary_path);

    Ok(written)
}

fn write_annotation_fasta<I>(path: &Path, records: I) -> std::io::Result<()>
where
    I: IntoIterator<Item = (String, String, String)>,
{
    let mut file = std::fs::File::create(path)?;
    for (identifier, description, sequence) in records {
        writeln!(file, ">{identifier} {description}")?;
        for chunk in sequence.as_bytes().chunks(80) {
            file.write_all(chunk)?;
            file.write_all(b"\n")?;
        }
    }
    Ok(())
}

fn format_record_description(record: &AnnotationRecord) -> String {
    let gene = &record.transcript.gene;
    format!(
        "gene:{} transcript:{} seqid:{} start:{} end:{} strand:{} class:{} family:{} complete:{}",
        gene.gene_id,
        record.transcript.transcript_id,
        record.transcript.seqid,
        gene.start0 + 1,
        gene.end0,
        record.transcript.strand.symbol(),
        record.arch.class,
        record.arch.family,
        record.arch.complete
    )
}

fn write_records_tsv(path: &Path, records: &[AnnotationRecord]) -> std::io::Result<()> {
    let refs: Vec<&AnnotationRecord> = records.iter().collect();
    write_record_refs_tsv(path, &refs)
}

fn write_record_refs_tsv(path: &Path, records: &[&AnnotationRecord]) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "#gene_id\tgene_name\tgene_biotype\ttranscript_id\ttranscript_name\tseqid\tgene_start\tgene_end\tstrand\tcds_count\tcds_length_bp\tprotein_length_aa\tclass\tfamily\tnterm\thas_nbarc\thas_tir\thas_lrr_any\thas_cterm_lrr\tcomplete\tn_domains\tdomains\tinternal_stop_count\thas_internal_stop\ttrailing_partial_bases\tphase_consistent"
    )?;
    for record in records {
        let gene = &record.transcript.gene;
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            gene.gene_id,
            gene.gene_name.as_deref().unwrap_or(""),
            gene.gene_biotype.as_deref().unwrap_or(""),
            record.transcript.transcript_id,
            record.transcript.transcript_name.as_deref().unwrap_or(""),
            record.transcript.seqid,
            gene.start0 + 1,
            gene.end0,
            record.transcript.strand.symbol(),
            record.transcript.cds.len(),
            record.transcript.cds_length_bp,
            record.protein_length(),
            record.arch.class,
            record.arch.family,
            record.arch.nterm,
            record.arch.has_nbarc,
            record.arch.has_tir,
            record.arch.has_lrr_any,
            record.arch.has_cterm_lrr,
            record.arch.complete,
            record.passed_domain_count(),
            record.domain_layout(),
            record.internal_stop_count,
            record.has_internal_stop,
            record.trailing_partial_bases,
            record.transcript.phase_consistent
        )?;
    }
    Ok(())
}

fn write_domains_tsv(path: &Path, records: &[AnnotationRecord]) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "#gene_id\ttranscript_id\tseqid\tstrand\tsource\tmodel\tmodel_name\tcategory\tprotein_start\tprotein_end\tcds_start\tcds_end\tgenomic_start\tgenomic_end\tgenomic_segments\tbitscore\tpvalue\tpassed"
    )?;
    for record in records {
        for domain in &record.domains {
            let pieces = domain
                .genomic_pieces
                .iter()
                .map(|piece| format!("{}:{}-{}", record.transcript.seqid, piece.start, piece.end))
                .collect::<Vec<_>>()
                .join(";");
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                record.transcript.gene.gene_id,
                record.transcript.transcript_id,
                record.transcript.seqid,
                record.transcript.strand.symbol(),
                domain.source,
                domain.model,
                domain.model_name,
                domain.category,
                domain.protein_start,
                domain.protein_end,
                domain.cds_start,
                domain.cds_end,
                domain.genomic_start,
                domain.genomic_end,
                pieces,
                domain
                    .bitscore
                    .map(|value| format!("{value:.1}"))
                    .unwrap_or_else(|| ".".to_string()),
                domain
                    .pvalue
                    .map(format_double_java)
                    .unwrap_or_else(|| ".".to_string()),
                domain.passed
            )?;
        }
    }
    Ok(())
}

fn write_skipped_tsv(path: &Path, skipped: &[SkippedAnnotation]) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(file, "#gene_id\ttranscript_id\tseqid\treason")?;
    for item in skipped {
        writeln!(
            file,
            "{}\t{}\t{}\t{}",
            item.gene_id, item.transcript_id, item.seqid, item.reason
        )?;
    }
    Ok(())
}

fn write_nlr_gff3(path: &Path, records: &[&AnnotationRecord]) -> std::io::Result<()> {
    let mut file = std::fs::File::create(path)?;
    writeln!(file, "##gff-version 3")?;
    writeln!(file, "##source FastNLR")?;
    for record in records {
        let gene = &record.transcript.gene;
        let mut attributes = format!(
            "ID={};Transcript={};Class={};Family={};Nterm={};DomainCount={};HasCtermLRR={};Complete={}",
            escape_gff3(&gene.gene_id),
            escape_gff3(&record.transcript.transcript_id),
            escape_gff3(&record.arch.class),
            escape_gff3(&record.arch.family),
            escape_gff3(&record.arch.nterm),
            record.passed_domain_count(),
            record.arch.has_cterm_lrr,
            record.arch.complete
        );
        if let Some(name) = &gene.gene_name {
            attributes.push_str(&format!(";Name={}", escape_gff3(name)));
        }
        if let Some(biotype) = &gene.gene_biotype {
            attributes.push_str(&format!(";gene_biotype={}", escape_gff3(biotype)));
        }
        writeln!(
            file,
            "{}\tFastNLR\tNLR_gene\t{}\t{}\t.\t{}\t.\t{}",
            record.transcript.seqid,
            gene.start0 + 1,
            gene.end0,
            record.transcript.strand.symbol(),
            attributes
        )?;
    }
    Ok(())
}

fn escape_gff3(value: &str) -> String {
    let mut output = String::with_capacity(value.len());
    for byte in value.bytes() {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'.' | b':' | b'_' | b'-' => {
                output.push(byte as char)
            }
            b' ' => output.push_str("%20"),
            _ => output.push_str(&format!("%{byte:02X}")),
        }
    }
    output
}

#[derive(Serialize)]
struct SummaryJson {
    mode: String,
    genome_file: String,
    gff_file: String,
    threads: usize,
    elapsed_s: f64,
    model_library_version: String,
    model_count: usize,
    domain_models: Vec<String>,
    user_hmm_files: Vec<String>,
    genes_with_cds: usize,
    genes_scanned: usize,
    skipped_genes: usize,
    nlr_genes: usize,
    complete_nlr_genes: usize,
    type_counts: BTreeMap<String, usize>,
    complete_type_counts: BTreeMap<String, usize>,
    family_counts: BTreeMap<String, usize>,
    skipped_reason_counts: BTreeMap<String, usize>,
}

fn summary_json(result: &AnnotationResult, meta: &AnnotationRunMeta) -> SummaryJson {
    let mut type_counts = BTreeMap::new();
    let mut complete_type_counts = BTreeMap::new();
    let mut family_counts = BTreeMap::new();
    for record in &result.records {
        if !record.arch.is_nlr {
            continue;
        }
        *type_counts.entry(record.arch.class.clone()).or_default() += 1;
        if record.arch.complete {
            *complete_type_counts
                .entry(record.arch.class.clone())
                .or_default() += 1;
        }
        *family_counts.entry(record.arch.family.clone()).or_default() += 1;
    }
    let mut skipped_reason_counts = BTreeMap::new();
    for item in &result.skipped {
        *skipped_reason_counts
            .entry(item.reason.clone())
            .or_default() += 1;
    }
    SummaryJson {
        mode: "annotation".to_string(),
        genome_file: meta.genome_name.clone(),
        gff_file: meta.gff_name.clone(),
        threads: meta.threads,
        elapsed_s: meta.elapsed.as_secs_f64(),
        model_library_version: meta.model_library_version.clone(),
        model_count: meta.model_count,
        domain_models: result.model_labels.clone(),
        user_hmm_files: meta.hmm_files.clone(),
        genes_with_cds: result.genes_with_cds,
        genes_scanned: result.records.len(),
        skipped_genes: result.skipped.len(),
        nlr_genes: result
            .records
            .iter()
            .filter(|record| record.arch.is_nlr)
            .count(),
        complete_nlr_genes: result
            .records
            .iter()
            .filter(|record| record.arch.is_nlr && record.arch.complete)
            .count(),
        type_counts,
        complete_type_counts,
        family_counts,
        skipped_reason_counts,
    }
}

fn write_summary_json(
    path: &Path,
    result: &AnnotationResult,
    meta: &AnnotationRunMeta,
) -> std::io::Result<()> {
    let report = summary_json(result, meta);
    let bytes = serde_json::to_vec_pretty(&report)
        .map_err(|error| std::io::Error::other(error.to_string()))?;
    std::fs::write(path, bytes)
}

fn write_summary_text(
    path: &Path,
    written: &[PathBuf],
    result: &AnnotationResult,
    meta: &AnnotationRunMeta,
) -> std::io::Result<()> {
    let report = summary_json(result, meta);
    let mut text = String::new();
    writeln!(text, "# FastNLR GFF3 annotation summary").unwrap();
    writeln!(text, "# input genome\t{}", meta.genome_name).unwrap();
    writeln!(text, "# input GFF3\t{}", meta.gff_name).unwrap();
    writeln!(text, "# threads\t{}", meta.threads).unwrap();
    writeln!(text, "# elapsed\t{:.2}s", meta.elapsed.as_secs_f64()).unwrap();
    writeln!(
        text,
        "# HMM model library\t{} ({} models)",
        meta.model_library_version, meta.model_count
    )
    .unwrap();
    writeln!(text, "# domain models\t{}", result.model_labels.join(", ")).unwrap();
    writeln!(text, "# genes with CDS\t{}", result.genes_with_cds).unwrap();
    writeln!(text, "# proteins scanned\t{}", result.records.len()).unwrap();
    writeln!(text, "# skipped genes\t{}", result.skipped.len()).unwrap();
    writeln!(text, "# NLR candidates\t{}", report.nlr_genes).unwrap();
    writeln!(
        text,
        "# Complete NLR candidates\t{}",
        report.complete_nlr_genes
    )
    .unwrap();
    writeln!(text, "# NLR types (count desc):").unwrap();
    write_sorted_counts(&mut text, &report.type_counts);
    writeln!(text, "# Complete NLR types (count desc):").unwrap();
    write_sorted_counts(&mut text, &report.complete_type_counts);
    writeln!(text, "# NLR families (count desc):").unwrap();
    write_sorted_counts(&mut text, &report.family_counts);
    writeln!(text, "# output files:").unwrap();
    for output in written {
        writeln!(text, "# {}", output.display()).unwrap();
    }
    std::fs::write(path, text)
}

fn write_sorted_counts(writer: &mut String, counts: &BTreeMap<String, usize>) {
    let mut values: Vec<(&String, &usize)> = counts.iter().collect();
    values.sort_by(|a, b| b.1.cmp(a.1).then_with(|| a.0.cmp(b.0)));
    for (name, count) in values {
        writeln!(writer, "# {count}\t{name}").unwrap();
    }
}

/// 供服务层和 GUI 生成一致的 GFF 模式元信息。
pub fn annotation_meta(
    config: &AnnotationConfig,
    elapsed: std::time::Duration,
    model_count: usize,
) -> AnnotationRunMeta {
    AnnotationRunMeta {
        genome_name: config
            .genome_fasta
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("genome")
            .to_string(),
        gff_name: config
            .gff3
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("annotation.gff3")
            .to_string(),
        threads: config.threads,
        elapsed,
        model_library_version: nlr_domain::MODEL_LIBRARY_VERSION.to_string(),
        model_count,
        hmm_files: config
            .hmm_files
            .iter()
            .map(|path| path.display().to_string())
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nlr_seq::gff3::{CdsSegment, GeneAnnotation};

    #[test]
    fn translation_trims_terminal_stop_and_replaces_internal_stop() {
        let translated = translate_cds("ATGTAAATG", 0);
        assert_eq!(translated.protein, "MXM");
        assert_eq!(translated.internal_stop_count, 1);
        assert_eq!(translated.trailing_partial_bases, 0);
    }

    #[test]
    fn translation_records_trailing_partial_codon() {
        let translated = translate_cds("ATGAA", 0);
        assert_eq!(translated.protein, "M");
        assert_eq!(translated.trailing_partial_bases, 2);
    }

    #[test]
    fn translation_skips_leading_partial_codon() {
        let translated = translate_cds("AATGAA", 1);
        assert_eq!(translated.protein, "M");
        assert_eq!(translated.trailing_partial_bases, 2);
    }

    #[test]
    fn maps_cross_exon_protein_interval_to_genomic_pieces() {
        let transcript = TranscriptModel {
            gene: GeneAnnotation {
                gene_id: "g1".to_string(),
                gene_name: None,
                gene_biotype: None,
                seqid: "chr1".to_string(),
                start0: 0,
                end0: 20,
                strand: Strand::Plus,
            },
            transcript_id: "t1".to_string(),
            transcript_name: None,
            seqid: "chr1".to_string(),
            start0: 0,
            end0: 20,
            strand: Strand::Plus,
            cds: vec![
                CdsSegment {
                    seqid: "chr1".to_string(),
                    start0: 0,
                    end0: 3,
                    phase: Some(0),
                },
                CdsSegment {
                    seqid: "chr1".to_string(),
                    start0: 10,
                    end0: 13,
                    phase: Some(0),
                },
            ],
            cds_length_bp: 6,
            phase_consistent: true,
        };
        let mapped = map_protein_interval(&transcript, 1, 2).unwrap();
        assert_eq!(mapped.0, 1);
        assert_eq!(mapped.1, 6);
        assert_eq!(
            mapped.4,
            vec![
                GenomicPiece { start: 1, end: 3 },
                GenomicPiece { start: 11, end: 13 }
            ]
        );
    }

    #[test]
    fn maps_protein_interval_with_leading_phase_offset() {
        let transcript = TranscriptModel {
            gene: GeneAnnotation {
                gene_id: "g1".to_string(),
                gene_name: None,
                gene_biotype: None,
                seqid: "chr1".to_string(),
                start0: 0,
                end0: 6,
                strand: Strand::Plus,
            },
            transcript_id: "t1".to_string(),
            transcript_name: None,
            seqid: "chr1".to_string(),
            start0: 0,
            end0: 6,
            strand: Strand::Plus,
            cds: vec![CdsSegment {
                seqid: "chr1".to_string(),
                start0: 0,
                end0: 6,
                phase: Some(1),
            }],
            cds_length_bp: 6,
            phase_consistent: true,
        };
        let mapped = map_protein_interval(&transcript, 1, 1).unwrap();
        assert_eq!((mapped.0, mapped.1), (2, 4));
        assert_eq!(mapped.4, vec![GenomicPiece { start: 2, end: 4 }]);
    }
}
