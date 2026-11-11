//! FastNLR 桌面端后端。
//!
//! Tauri 命令只负责参数校验、作业生命周期和结果适配；实际注释流水线全部复用
//! `nlr-service` / `nlr-cli`，不在桌面端复制扫描或组装逻辑。
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use std::collections::BTreeMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use nlr_report::SummaryReport;
use nlr_service::{AnnotationJob, GenomeJob, JobSnapshot, RunConfig, Service};
use serde::{Deserialize, Serialize};
use tauri::State;

#[derive(Clone)]
struct AppState {
    service: Service,
}

/// 桌面端运行环境的只读信息。
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct RuntimeInfo {
    version: String,
    os: String,
    arch: String,
    logical_cpus: usize,
}

/// 输入模式与结果模式共用同一枚举，前端只需维护一个 discriminant。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
enum InputMode {
    Genome,
    Annotation,
}

/// 通过文件选择器选中后的输入文件信息。
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct InputInfo {
    path: String,
    name: String,
    prefix: String,
    size_bytes: u64,
    compressed: bool,
}

/// 后端返回给前端的单个输出条目。
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct OutputFile {
    name: String,
    path: String,
    size_bytes: u64,
    is_dir: bool,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct GenomeJobRequest {
    input_path: String,
    output_dir: String,
    #[serde(default)]
    threads: Option<usize>,
    #[serde(default)]
    flank: Option<u64>,
    #[serde(default)]
    seqs_per_thread: Option<usize>,
    #[serde(default)]
    checkpoint_dir: Option<String>,
    #[serde(default)]
    motif_accept_p: Option<f64>,
    #[serde(default)]
    motif_prelim_p: Option<f64>,
    #[serde(default)]
    tir_accept_p: Option<f64>,
    #[serde(default)]
    flexible_seed: Option<usize>,
    #[serde(default)]
    no_flexible_seed: bool,
    #[serde(default)]
    require_ploop: bool,
    #[serde(default)]
    relaxed_seed: Option<usize>,
    #[serde(default)]
    keep_partial_loci: bool,
    #[serde(default)]
    mot_file: Option<String>,
    #[serde(default)]
    store_file: Option<String>,
    #[serde(default)]
    motif_categories: Vec<String>,
    #[serde(default)]
    extra_seeds: Vec<String>,
    #[serde(default)]
    extra_signatures: Vec<String>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct AnnotationJobRequest {
    input_path: String,
    gff_path: String,
    output_dir: String,
    #[serde(default)]
    threads: Option<usize>,
    #[serde(default)]
    batch: Option<usize>,
    // Retained only to reject pre-1.9.1 PWM fields explicitly.
    #[serde(default)]
    motif_accept_p: Option<f64>,
    #[serde(default)]
    motif_prelim_p: Option<f64>,
    #[serde(default)]
    hmm_files: Vec<String>,
    #[serde(default)]
    mot_file: Option<String>,
    #[serde(default)]
    store_file: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct ResultRequest {
    output_dir: String,
    prefix: String,
    mode: InputMode,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct AnnotationRecordView {
    gene_id: String,
    gene_name: String,
    gene_biotype: String,
    transcript_id: String,
    transcript_name: String,
    seqid: String,
    gene_start: u64,
    gene_end: u64,
    strand: String,
    cds_count: usize,
    cds_length_bp: u64,
    protein_length: usize,
    class: String,
    nterm: String,
    has_nbarc: bool,
    has_tir: bool,
    has_lrr_any: bool,
    has_cterm_lrr: bool,
    complete: bool,
    is_nlr: bool,
    domain_count: usize,
    domains: String,
    family: String,
    internal_stop_count: usize,
    has_internal_stop: bool,
    trailing_partial_bases: u8,
    phase_consistent: bool,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct AnnotationDomainView {
    gene_id: String,
    transcript_id: String,
    seqid: String,
    strand: String,
    source: String,
    model: String,
    model_name: String,
    category: String,
    protein_start: usize,
    protein_end: usize,
    cds_start: usize,
    cds_end: usize,
    genomic_start: u64,
    genomic_end: u64,
    genomic_segments: String,
    bitscore: Option<f32>,
    pvalue: Option<f64>,
    passed: bool,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct AnnotationResults {
    proteins_scanned: usize,
    nlr_total: usize,
    complete_nlr_total: usize,
    genes_with_cds: usize,
    skipped_genes: usize,
    elapsed: Option<String>,
    model_labels: Vec<String>,
    model_library_version: String,
    model_count: usize,
    type_counts: BTreeMap<String, usize>,
    complete_type_counts: BTreeMap<String, usize>,
    family_counts: BTreeMap<String, usize>,
    records: Vec<AnnotationRecordView>,
    domains: Vec<AnnotationDomainView>,
}

#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
struct RunResults {
    mode: InputMode,
    prefix: String,
    output_dir: String,
    outputs: Vec<OutputFile>,
    genome: Option<SummaryReport>,
    annotation: Option<AnnotationResults>,
}

/// 返回桌面端版本、目标和可用逻辑 CPU，供设置页与默认线程数使用。
#[tauri::command]
fn get_runtime_info() -> RuntimeInfo {
    RuntimeInfo {
        version: env!("CARGO_PKG_VERSION").to_string(),
        os: std::env::consts::OS.to_string(),
        arch: std::env::consts::ARCH.to_string(),
        logical_cpus: std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1),
    }
}

/// 检查输入文件并推导与 CLI 一致的文件名前缀。
#[tauri::command]
fn inspect_input(path: String, mode: InputMode) -> Result<InputInfo, String> {
    let input = PathBuf::from(&path);
    validate_input(&input, mode)?;
    let metadata = fs::metadata(&input).map_err(|e| format!("cannot read input metadata: {e}"))?;
    let name = input
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| "input filename is not valid UTF-8".to_string())?
        .to_string();
    Ok(InputInfo {
        path: path.clone(),
        prefix: derive_prefix(&path),
        name,
        size_bytes: metadata.len(),
        compressed: path.to_ascii_lowercase().ends_with(".gz"),
    })
}

/// 提交基因组注释作业。参数构造和校验与 CLI 使用同一 `RunConfig`。
#[tauri::command]
fn start_genome_job(state: State<'_, AppState>, request: GenomeJobRequest) -> Result<u64, String> {
    let job = build_genome_job(request)?;
    Ok(state.service.submit_genome(job))
}

/// 提交 GFF3 驱动蛋白注释作业。
#[tauri::command]
fn start_annotation_job(
    state: State<'_, AppState>,
    request: AnnotationJobRequest,
) -> Result<u64, String> {
    let job = build_annotation_job(request)?;
    Ok(state.service.submit_annotation(job))
}

/// 读取作业快照；前端用约 500 ms 的间隔轮询，避免引入额外事件桥接复杂度。
#[tauri::command]
fn get_job(state: State<'_, AppState>, job_id: u64) -> Result<JobSnapshot, String> {
    state
        .service
        .snapshot(job_id)
        .ok_or_else(|| format!("job {job_id} does not exist"))
}

#[tauri::command]
fn cancel_job(state: State<'_, AppState>, job_id: u64) -> bool {
    state.service.cancel(job_id)
}

/// 判断同一前缀是否已有完成结果，前端据此要求用户显式确认覆盖。
#[tauri::command]
fn has_existing_results(output_dir: String, prefix: String, mode: InputMode) -> bool {
    let output_dir = Path::new(&output_dir);
    let filename = match mode {
        InputMode::Genome => format!("{prefix}.summary.json"),
        InputMode::Annotation => format!("{prefix}.gff.summary.json"),
    };
    output_dir.join(filename).is_file()
}

/// 从现有输出目录加载结果，支持应用重启后的历史记录回看。
#[tauri::command]
async fn load_results(request: ResultRequest) -> Result<RunResults, String> {
    tauri::async_runtime::spawn_blocking(move || load_results_blocking(&request))
        .await
        .map_err(|e| format!("result loader task failed: {e}"))?
}

fn build_genome_job(request: GenomeJobRequest) -> Result<GenomeJob, String> {
    let input_path = PathBuf::from(&request.input_path);
    validate_input(&input_path, InputMode::Genome)?;
    if request.relaxed_seed.is_some() && !request.no_flexible_seed {
        return Err(
            "relaxed seeding requires disabling flexible seeding (noFlexibleSeed)".to_string(),
        );
    }

    let mut config = RunConfig::new(
        input_path.clone(),
        optional_path(request.mot_file),
        optional_path(request.store_file),
    );
    if let Some(threads) = request.threads {
        config.threads = threads;
    }
    if let Some(batch) = request.seqs_per_thread {
        config.seqs_per_thread = batch;
    }
    if let Some(value) = request.motif_accept_p {
        config.motif_accept_p = value;
    }
    if let Some(value) = request.motif_prelim_p {
        config.motif_prelim_p = value;
    }
    if let Some(value) = request.tir_accept_p {
        config.tir_accept_p = value;
    }
    config.checkpoint_dir = optional_path(request.checkpoint_dir);
    config.assemble.relaxed_seed_min_nbarc = request.relaxed_seed;
    config.flexible_seed = nlr_cli::resolve_flexible_seed(
        request.flexible_seed,
        request.no_flexible_seed,
        request.require_ploop,
    );
    config.locus_filter = if request.keep_partial_loci {
        None
    } else {
        Some(nlr_cli::LocusFilterConfig::high_precision_default())
    };
    config.motif_categories = request.motif_categories;
    config.extra_seeds = request.extra_seeds;
    config.extra_signatures = request.extra_signatures;
    config
        .validate()
        .map_err(|e| format!("invalid genome parameters: {e}"))?;
    let output_dir = PathBuf::from(&request.output_dir);
    fs::create_dir_all(&output_dir).map_err(|e| format!("cannot create output directory: {e}"))?;

    Ok(GenomeJob {
        prefix: derive_prefix(&request.input_path),
        config,
        output_dir,
        flank: request.flank.unwrap_or(2000),
    })
}

fn build_annotation_job(request: AnnotationJobRequest) -> Result<AnnotationJob, String> {
    if request.motif_accept_p.is_some()
        || request.motif_prelim_p.is_some()
        || request.mot_file.is_some()
        || request.store_file.is_some()
    {
        return Err(
            "GFF3 annotation mode uses HMM-only classification; PWM fields are not supported"
                .to_string(),
        );
    }
    let input_path = PathBuf::from(&request.input_path);
    validate_input(&input_path, InputMode::Annotation)?;
    let gff_path = PathBuf::from(&request.gff_path);
    validate_gff(&gff_path)?;
    let mut config = nlr_cli::AnnotationConfig::new(input_path, gff_path);
    if let Some(threads) = request.threads {
        config.threads = threads;
    }
    if let Some(batch) = request.batch {
        config.batch = batch;
    }
    config.hmm_files = request.hmm_files.into_iter().map(PathBuf::from).collect();
    config
        .validate()
        .map_err(|e| format!("invalid GFF3 annotation parameters: {e}"))?;
    let output_dir = PathBuf::from(&request.output_dir);
    fs::create_dir_all(&output_dir).map_err(|e| format!("cannot create output directory: {e}"))?;

    Ok(AnnotationJob {
        prefix: derive_prefix(&request.input_path),
        config,
        output_dir,
    })
}

fn load_results_blocking(request: &ResultRequest) -> Result<RunResults, String> {
    if request.prefix.is_empty()
        || request.prefix.contains('/')
        || request.prefix.contains('\\')
        || request.prefix == "."
        || request.prefix == ".."
    {
        return Err("invalid output prefix".to_string());
    }
    let output_dir = PathBuf::from(&request.output_dir);
    if !output_dir.is_dir() {
        return Err(format!(
            "output directory does not exist: {}",
            output_dir.display()
        ));
    }

    let (genome, annotation) = match request.mode {
        InputMode::Genome => (
            Some(read_genome_results(&output_dir, &request.prefix)?),
            None,
        ),
        InputMode::Annotation => (
            None,
            Some(read_annotation_results(&output_dir, &request.prefix)?),
        ),
    };

    Ok(RunResults {
        mode: request.mode,
        prefix: request.prefix.clone(),
        output_dir: request.output_dir.clone(),
        outputs: list_outputs(&output_dir, &request.prefix)?,
        genome,
        annotation,
    })
}

fn read_genome_results(output_dir: &Path, prefix: &str) -> Result<SummaryReport, String> {
    let path = output_dir.join(format!("{prefix}.summary.json"));
    let file = fs::File::open(&path)
        .map_err(|e| format!("cannot read genome summary {}: {e}", path.display()))?;
    serde_json::from_reader(BufReader::new(file))
        .map_err(|e| format!("cannot parse genome summary {}: {e}", path.display()))
}

#[derive(Debug, Deserialize)]
struct AnnotationSummaryFile {
    genes_with_cds: usize,
    genes_scanned: usize,
    skipped_genes: usize,
    nlr_genes: usize,
    #[serde(default)]
    complete_nlr_genes: usize,
    elapsed_s: f64,
    domain_models: Vec<String>,
    #[serde(default)]
    model_library_version: String,
    #[serde(default)]
    model_count: usize,
    type_counts: BTreeMap<String, usize>,
    #[serde(default)]
    complete_type_counts: BTreeMap<String, usize>,
    family_counts: BTreeMap<String, usize>,
}

fn read_annotation_results(output_dir: &Path, prefix: &str) -> Result<AnnotationResults, String> {
    let genes_path = output_dir.join(format!("{prefix}.gff.genes.tsv"));
    let domain_path = output_dir.join(format!("{prefix}.gff.domains.tsv"));
    let summary_path = output_dir.join(format!("{prefix}.gff.summary.json"));

    let mut records = parse_annotation_records(&genes_path)?;
    records.retain(|record| record.is_nlr);
    let domains = parse_annotation_domains(&domain_path)?;
    let summary_file = fs::File::open(&summary_path).map_err(|e| {
        format!(
            "cannot read annotation summary {}: {e}",
            summary_path.display()
        )
    })?;
    let summary: AnnotationSummaryFile = serde_json::from_reader(BufReader::new(summary_file))
        .map_err(|e| {
            format!(
                "cannot parse annotation summary {}: {e}",
                summary_path.display()
            )
        })?;

    Ok(AnnotationResults {
        proteins_scanned: summary.genes_scanned,
        nlr_total: summary.nlr_genes,
        complete_nlr_total: summary.complete_nlr_genes,
        genes_with_cds: summary.genes_with_cds,
        skipped_genes: summary.skipped_genes,
        elapsed: Some(format!("{:.2}s", summary.elapsed_s)),
        model_labels: summary.domain_models,
        model_library_version: summary.model_library_version,
        model_count: summary.model_count,
        type_counts: summary.type_counts,
        complete_type_counts: summary.complete_type_counts,
        family_counts: summary.family_counts,
        records,
        domains,
    })
}

fn parse_annotation_records(path: &Path) -> Result<Vec<AnnotationRecordView>, String> {
    let file = fs::File::open(path)
        .map_err(|e| format!("cannot read annotation records {}: {e}", path.display()))?;
    let mut records = Vec::new();
    for (line_number, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if !matches!(columns.len(), 24 | 26) {
            return Err(format!(
                "{}:{} expected 24 or 26 columns, found {}",
                path.display(),
                line_number + 1,
                columns.len()
            ));
        }
        let has_nbarc = parse_bool(columns[15]);
        let has_tir = parse_bool(columns[16]);
        let has_lrr_any = parse_bool(columns[17]);
        let has_cterm_lrr = if columns.len() == 26 {
            parse_bool(columns[18])
        } else {
            false
        };
        let complete = if columns.len() == 26 {
            parse_bool(columns[19])
        } else {
            false
        };
        records.push(AnnotationRecordView {
            gene_id: columns[0].to_string(),
            gene_name: columns[1].to_string(),
            gene_biotype: columns[2].to_string(),
            transcript_id: columns[3].to_string(),
            transcript_name: columns[4].to_string(),
            seqid: columns[5].to_string(),
            gene_start: parse_u64(columns[6], path, line_number)?,
            gene_end: parse_u64(columns[7], path, line_number)?,
            strand: columns[8].to_string(),
            cds_count: parse_usize(columns[9], path, line_number)?,
            cds_length_bp: parse_u64(columns[10], path, line_number)?,
            protein_length: parse_usize(columns[11], path, line_number)?,
            class: columns[12].to_string(),
            family: columns[13].to_string(),
            nterm: columns[14].to_string(),
            has_nbarc,
            has_tir,
            has_lrr_any,
            has_cterm_lrr,
            complete,
            is_nlr: has_nbarc || (has_tir && has_cterm_lrr),
            domain_count: parse_usize(
                columns[if columns.len() == 26 { 20 } else { 18 }],
                path,
                line_number,
            )?,
            domains: columns[if columns.len() == 26 { 21 } else { 19 }].to_string(),
            internal_stop_count: parse_usize(
                columns[if columns.len() == 26 { 22 } else { 20 }],
                path,
                line_number,
            )?,
            has_internal_stop: parse_bool(columns[if columns.len() == 26 { 23 } else { 21 }]),
            trailing_partial_bases: columns[if columns.len() == 26 { 24 } else { 22 }]
                .parse()
                .map_err(|_| {
                    format!(
                        "{}:{} invalid integer '{}'",
                        path.display(),
                        line_number + 1,
                        columns[if columns.len() == 26 { 24 } else { 22 }]
                    )
                })?,
            phase_consistent: parse_bool(columns[if columns.len() == 26 { 25 } else { 23 }]),
        });
    }
    Ok(records)
}

fn parse_annotation_domains(path: &Path) -> Result<Vec<AnnotationDomainView>, String> {
    let file = fs::File::open(path)
        .map_err(|e| format!("cannot read annotation domains {}: {e}", path.display()))?;
    let mut domains = Vec::new();
    for (line_number, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|e| format!("cannot read {}: {e}", path.display()))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() != 18 {
            return Err(format!(
                "{}:{} expected 18 columns, found {}",
                path.display(),
                line_number + 1,
                columns.len()
            ));
        }
        domains.push(AnnotationDomainView {
            gene_id: columns[0].to_string(),
            transcript_id: columns[1].to_string(),
            seqid: columns[2].to_string(),
            strand: columns[3].to_string(),
            source: columns[4].to_string(),
            model: columns[5].to_string(),
            model_name: columns[6].to_string(),
            category: columns[7].to_string(),
            protein_start: parse_usize(columns[8], path, line_number)?,
            protein_end: parse_usize(columns[9], path, line_number)?,
            cds_start: parse_usize(columns[10], path, line_number)?,
            cds_end: parse_usize(columns[11], path, line_number)?,
            genomic_start: parse_u64(columns[12], path, line_number)?,
            genomic_end: parse_u64(columns[13], path, line_number)?,
            genomic_segments: columns[14].to_string(),
            bitscore: columns[15].parse().ok(),
            pvalue: columns[16].parse().ok(),
            passed: parse_bool(columns[17]),
        });
    }
    Ok(domains)
}

fn list_outputs(output_dir: &Path, prefix: &str) -> Result<Vec<OutputFile>, String> {
    let mut outputs = Vec::new();
    let entries = fs::read_dir(output_dir)
        .map_err(|e| format!("cannot list output directory {}: {e}", output_dir.display()))?;
    for entry in entries {
        let entry = entry.map_err(|e| format!("cannot read output entry: {e}"))?;
        let name = entry.file_name().to_string_lossy().to_string();
        if !name.starts_with(prefix) {
            continue;
        }
        let path = entry.path();
        let metadata = entry
            .metadata()
            .map_err(|e| format!("cannot read metadata for {}: {e}", path.display()))?;
        outputs.push(OutputFile {
            name,
            path: path.display().to_string(),
            size_bytes: if metadata.is_file() {
                metadata.len()
            } else {
                0
            },
            is_dir: metadata.is_dir(),
        });
    }
    outputs.sort_by(|left, right| left.name.cmp(&right.name));
    Ok(outputs)
}

fn validate_input(path: &Path, _mode: InputMode) -> Result<(), String> {
    if !path.is_file() {
        return Err(format!("input file does not exist: {}", path.display()));
    }
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| "input filename is not valid UTF-8".to_string())?
        .to_ascii_lowercase();
    let uncompressed = name.strip_suffix(".gz").unwrap_or(&name);
    let extension = Path::new(uncompressed)
        .extension()
        .and_then(|extension| extension.to_str())
        .unwrap_or_default();
    if !matches!(extension, "fa" | "fasta" | "fna" | "fas" | "ffn") {
        return Err(format!(
            "unsupported genome FASTA extension '.{extension}'; expected .fa, .fasta, .fna, .fas, .ffn, or .gz"
        ));
    }
    Ok(())
}

fn validate_gff(path: &Path) -> Result<(), String> {
    if !path.is_file() {
        return Err(format!("GFF3 file does not exist: {}", path.display()));
    }
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| "GFF3 filename is not valid UTF-8".to_string())?
        .to_ascii_lowercase();
    let uncompressed = name.strip_suffix(".gz").unwrap_or(&name);
    let extension = Path::new(uncompressed)
        .extension()
        .and_then(|extension| extension.to_str())
        .unwrap_or_default();
    if !matches!(extension, "gff" | "gff3") {
        return Err(format!(
            "unsupported GFF3 extension '.{extension}'; expected .gff, .gff3, or .gz"
        ));
    }
    Ok(())
}

fn derive_prefix(input: &str) -> String {
    let filename = input
        .rsplit(['/', '\\'])
        .next()
        .filter(|name| !name.is_empty())
        .unwrap_or("genome");
    let uncompressed = if filename.to_ascii_lowercase().ends_with(".gz") {
        &filename[..filename.len() - 3]
    } else {
        filename
    };
    let path = Path::new(uncompressed);
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    if matches!(extension.as_str(), "fa" | "fasta" | "fna" | "fas" | "ffn") {
        path.file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or("genome")
            .to_string()
    } else {
        uncompressed.to_string()
    }
}

fn optional_path(value: Option<String>) -> Option<PathBuf> {
    value.and_then(|value| {
        let trimmed = value.trim();
        (!trimmed.is_empty()).then(|| PathBuf::from(trimmed))
    })
}

fn parse_usize(value: &str, path: &Path, line_number: usize) -> Result<usize, String> {
    value.parse().map_err(|_| {
        format!(
            "{}:{} invalid integer '{value}'",
            path.display(),
            line_number + 1
        )
    })
}

fn parse_u64(value: &str, path: &Path, line_number: usize) -> Result<u64, String> {
    value.parse().map_err(|_| {
        format!(
            "{}:{} invalid integer '{value}'",
            path.display(),
            line_number + 1
        )
    })
}

fn parse_bool(value: &str) -> bool {
    matches!(value.trim(), "true" | "TRUE" | "1")
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .plugin(tauri_plugin_store::Builder::new().build())
        .manage(AppState {
            service: Service::new(),
        })
        .invoke_handler(tauri::generate_handler![
            get_runtime_info,
            inspect_input,
            start_genome_job,
            start_annotation_job,
            get_job,
            cancel_job,
            has_existing_results,
            load_results
        ])
        .run(tauri::generate_context!())
        .expect("error while running FastNLR desktop app");
}

#[cfg(test)]
mod tests {
    use super::*;
    use nlr_service::AnnotationConfig;

    fn genome_request(input: &Path, output: &Path) -> GenomeJobRequest {
        GenomeJobRequest {
            input_path: input.display().to_string(),
            output_dir: output.display().to_string(),
            threads: Some(2),
            flank: Some(2000),
            seqs_per_thread: Some(100),
            checkpoint_dir: None,
            motif_accept_p: Some(0.00001),
            motif_prelim_p: Some(0.0001),
            tir_accept_p: Some(0.0000001),
            flexible_seed: Some(3),
            no_flexible_seed: false,
            require_ploop: false,
            relaxed_seed: None,
            keep_partial_loci: false,
            mot_file: None,
            store_file: None,
            motif_categories: Vec::new(),
            extra_seeds: Vec::new(),
            extra_signatures: Vec::new(),
        }
    }

    #[test]
    fn derives_prefix_for_plain_and_gzip_fasta() {
        assert_eq!(derive_prefix("/tmp/genome.fa"), "genome");
        assert_eq!(derive_prefix("/tmp/genome.fasta.gz"), "genome");
        assert_eq!(derive_prefix("C:\\data\\genome.ffn"), "genome");
    }

    #[test]
    fn rejects_missing_input_and_creates_genome_job() {
        let missing = genome_request(Path::new("/tmp/does-not-exist.fa"), Path::new("/tmp/out"));
        assert!(build_genome_job(missing).is_err());
    }

    #[test]
    fn validates_conflicting_seeding_modes_before_creating_output() {
        let tmp = tempfile::TempDir::new().unwrap();
        let input = tmp.path().join("genome.fa");
        fs::write(&input, ">chr1\nATGATGATG\n").unwrap();
        let output = tmp.path().join("results");
        let mut request = genome_request(&input, &output);
        request.relaxed_seed = Some(3);
        let error = build_genome_job(request).unwrap_err();
        assert!(error.contains("requires disabling flexible seeding"));
        assert!(!output.exists());
    }

    #[test]
    fn rejects_pwm_fields_in_annotation_mode() {
        let request = AnnotationJobRequest {
            input_path: "/tmp/genome.fa".to_string(),
            gff_path: "/tmp/annotation.gff3".to_string(),
            output_dir: "/tmp/fastnlr-annotation".to_string(),
            threads: Some(1),
            batch: Some(10),
            motif_accept_p: Some(0.001),
            motif_prelim_p: None,
            hmm_files: Vec::new(),
            mot_file: None,
            store_file: None,
        };
        let error = build_annotation_job(request).unwrap_err();
        assert!(error.contains("HMM-only"));
    }

    #[test]
    fn annotation_pipeline_output_is_loaded_into_desktop_results() {
        let tmp = tempfile::TempDir::new().unwrap();
        let genome = tmp.path().join("genome.fa");
        fs::write(&genome, ">chr1\nATGAAAACCGGGTTTAAATAA\n").unwrap();
        let gff = tmp.path().join("annotation.gff3");
        fs::write(
            &gff,
            concat!(
                "##gff-version 3\n",
                "chr1\ttest\tgene\t1\t21\t.\t+\t.\tID=g1;Name=test;gene_biotype=protein_coding\n",
                "chr1\ttest\tmRNA\t1\t21\t.\t+\t.\tID=t1;Parent=g1\n",
                "chr1\ttest\tCDS\t1\t21\t.\t+\t0\tParent=t1\n"
            ),
        )
        .unwrap();
        let output_dir = tmp.path().join("results");
        let mut config = AnnotationConfig::new(genome, gff);
        config.threads = 1;
        config.batch = 1;
        let job = AnnotationJob {
            config,
            output_dir: output_dir.clone(),
            prefix: "genome".to_string(),
        };
        nlr_service::run_annotation_blocking(&job).unwrap();

        let results = load_results_blocking(&ResultRequest {
            output_dir: output_dir.display().to_string(),
            prefix: "genome".to_string(),
            mode: InputMode::Annotation,
        })
        .unwrap();
        let annotation = results.annotation.unwrap();
        assert_eq!(annotation.proteins_scanned, 1);
        assert!(annotation.records.is_empty());
        assert!(results
            .outputs
            .iter()
            .any(|file| file.name.ends_with(".gff.summary.json")));
    }

    #[test]
    fn annotation_record_view_derives_nlr_candidate_status() {
        let tmp = tempfile::TempDir::new().unwrap();
        let path = tmp.path().join("genes.tsv");
        fs::write(
            &path,
            concat!(
                "#header\n",
                "g1\tname\tprotein_coding\tt1\t\tchr1\t1\t100\t+\t1\t90\t30\t",
                "CC-NBARC-LRR\tCNL\tCC\ttrue\tfalse\ttrue\ttrue\ttrue\t3\t",
                "PF00931:5-20\t0\tfalse\t0\ttrue\n",
                "g2\tname\tprotein_coding\tt2\t\tchr1\t200\t300\t+\t1\t90\t30\t",
                "Other\tOther\tNone\tfalse\tfalse\tfalse\tfalse\tfalse\t0\t\t",
                "0\tfalse\t0\ttrue\n"
            ),
        )
        .unwrap();
        let records = parse_annotation_records(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert!(records[0].is_nlr);
        assert!(!records[1].is_nlr);
    }
}
