//! nlr-service — 供桌面版 / Web 服务复用的作业执行层。
//!
//! 提供统一的作业生命周期：提交、查询进度、取消与结果收集。内部在独立线程中调用
//! `nlr_cli` 的 `run_tracked` / `run_protein_tracked`，不复制扫描与组装逻辑。
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use serde::Serialize;

pub use nlr_cli::{ProteinConfig, RunConfig};

/// 作业类型。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum JobKind {
    Genome,
    Protein,
}

/// 作业状态。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum JobStatus {
    Pending,
    Running,
    Completed,
    Failed,
    Canceled,
}

/// 作业执行阶段，用于 GUI 展示比 Running 更具体的工作状态。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum JobPhase {
    Pending,
    Scanning,
    Writing,
    Completed,
    Failed,
    Canceled,
}

/// 作业状态快照（可序列化给前端）。
#[derive(Debug, Clone, Serialize)]
pub struct JobSnapshot {
    pub id: u64,
    pub kind: JobKind,
    pub status: JobStatus,
    pub phase: JobPhase,
    pub scanned: u64,
    pub estimated_total: u64,
    pub error: Option<String>,
}

/// 基因组模式作业规格。
#[derive(Debug, Clone)]
pub struct GenomeJob {
    pub config: RunConfig,
    pub output_dir: PathBuf,
    pub prefix: String,
    pub flank: u64,
}

/// 蛋白模式作业规格。
#[derive(Debug, Clone)]
pub struct ProteinJob {
    pub config: ProteinConfig,
    pub output_dir: PathBuf,
    pub prefix: String,
}

struct JobState {
    kind: JobKind,
    status: Mutex<JobStatus>,
    phase: Mutex<JobPhase>,
    scanned: AtomicU64,
    estimated_total: AtomicU64,
    error: Mutex<Option<String>>,
    outputs: Mutex<Option<Vec<PathBuf>>>,
    cancel: AtomicBool,
}

/// 线程安全的作业注册表；每个作业在独立线程中执行。
#[derive(Clone)]
pub struct Service {
    jobs: Arc<Mutex<HashMap<u64, Arc<JobState>>>>,
    next_id: Arc<AtomicU64>,
}

impl Service {
    pub fn new() -> Self {
        Service {
            jobs: Arc::new(Mutex::new(HashMap::new())),
            next_id: Arc::new(AtomicU64::new(0)),
        }
    }

    /// 提交基因组作业，返回作业 id。
    pub fn submit_genome(&self, job: GenomeJob) -> u64 {
        let id = self.next_id.fetch_add(1, Ordering::SeqCst) + 1;
        let state = Arc::new(JobState {
            kind: JobKind::Genome,
            status: Mutex::new(JobStatus::Pending),
            phase: Mutex::new(JobPhase::Pending),
            scanned: AtomicU64::new(0),
            estimated_total: AtomicU64::new(estimate_total_fragments(
                &job.config.input_fasta,
                job.config.fragment_length,
                job.config.overlap,
            )),
            error: Mutex::new(None),
            outputs: Mutex::new(None),
            cancel: AtomicBool::new(false),
        });
        self.jobs.lock().unwrap().insert(id, state.clone());
        let service = self.clone();
        std::thread::spawn(move || service.run_genome_job(id, state, job));
        id
    }

    /// 提交蛋白作业，返回作业 id。
    pub fn submit_protein(&self, job: ProteinJob) -> u64 {
        let id = self.next_id.fetch_add(1, Ordering::SeqCst) + 1;
        let state = Arc::new(JobState {
            kind: JobKind::Protein,
            status: Mutex::new(JobStatus::Pending),
            phase: Mutex::new(JobPhase::Pending),
            scanned: AtomicU64::new(0),
            estimated_total: AtomicU64::new(0),
            error: Mutex::new(None),
            outputs: Mutex::new(None),
            cancel: AtomicBool::new(false),
        });
        self.jobs.lock().unwrap().insert(id, state.clone());
        let service = self.clone();
        std::thread::spawn(move || service.run_protein_job(id, state, job));
        id
    }

    /// 请求取消作业；返回作业是否存在。
    pub fn cancel(&self, id: u64) -> bool {
        let jobs = self.jobs.lock().unwrap();
        match jobs.get(&id) {
            Some(state) => {
                state.cancel.store(true, Ordering::SeqCst);
                true
            }
            None => false,
        }
    }

    /// 读取作业状态快照。
    pub fn snapshot(&self, id: u64) -> Option<JobSnapshot> {
        let state = {
            let jobs = self.jobs.lock().unwrap();
            jobs.get(&id)?.clone()
        };
        let status = *state.status.lock().unwrap();
        let phase = *state.phase.lock().unwrap();
        let error = state.error.lock().unwrap().clone();
        Some(JobSnapshot {
            id,
            kind: state.kind,
            status,
            phase,
            scanned: state.scanned.load(Ordering::SeqCst),
            estimated_total: state.estimated_total.load(Ordering::SeqCst),
            error,
        })
    }

    /// 读取作业的输出文件路径（完成前为 None）。
    pub fn outputs(&self, id: u64) -> Option<Vec<PathBuf>> {
        let state = {
            let jobs = self.jobs.lock().unwrap();
            jobs.get(&id)?.clone()
        };
        let outputs = state.outputs.lock().unwrap().clone();
        outputs
    }

    fn run_genome_job(&self, id: u64, state: Arc<JobState>, job: GenomeJob) {
        *state.status.lock().unwrap() = JobStatus::Running;
        *state.phase.lock().unwrap() = JobPhase::Scanning;
        let state_for_progress = state.clone();
        let state_for_phase = state.clone();
        let outcome = execute_genome(
            &job,
            move |fragments| {
                state_for_progress
                    .scanned
                    .store(fragments, Ordering::SeqCst);
                state_for_progress.cancel.load(Ordering::SeqCst)
            },
            move |phase| *state_for_phase.phase.lock().unwrap() = phase,
        );
        finalize(state, outcome);
        let _ = id;
    }

    fn run_protein_job(&self, id: u64, state: Arc<JobState>, job: ProteinJob) {
        *state.status.lock().unwrap() = JobStatus::Running;
        *state.phase.lock().unwrap() = JobPhase::Scanning;
        let state_for_progress = state.clone();
        let state_for_phase = state.clone();
        let outcome = execute_protein(
            &job,
            move |proteins| {
                state_for_progress.scanned.store(proteins, Ordering::SeqCst);
                state_for_progress.cancel.load(Ordering::SeqCst)
            },
            move |phase| *state_for_phase.phase.lock().unwrap() = phase,
        );
        finalize(state, outcome);
        let _ = id;
    }
}

impl Default for Service {
    fn default() -> Self {
        Self::new()
    }
}

/// 同步执行基因组作业（阻塞调用用）；返回输出文件路径。
pub fn run_genome_blocking(job: &GenomeJob) -> Result<Vec<PathBuf>, String> {
    match execute_genome(job, |_| false, |_| {}) {
        Ok(paths) => Ok(paths),
        Err(ExecutionError::Canceled) => Err("job canceled".to_string()),
        Err(ExecutionError::Failed(message)) => Err(message),
    }
}

/// 同步执行蛋白作业（阻塞调用用）；返回输出文件路径。
pub fn run_protein_blocking(job: &ProteinJob) -> Result<Vec<PathBuf>, String> {
    match execute_protein(job, |_| false, |_| {}) {
        Ok(paths) => Ok(paths),
        Err(ExecutionError::Canceled) => Err("job canceled".to_string()),
        Err(ExecutionError::Failed(message)) => Err(message),
    }
}

enum ExecutionError {
    Canceled,
    Failed(String),
}

fn execute_genome(
    job: &GenomeJob,
    mut progress: impl FnMut(u64) -> bool,
    mut phase: impl FnMut(JobPhase),
) -> Result<Vec<PathBuf>, ExecutionError> {
    let mut canceled = false;
    let mut result = nlr_cli::run_tracked(&job.config, |p| {
        if progress(p.scanned_fragments) {
            canceled = true;
            nlr_cli::RunControl::Break
        } else {
            nlr_cli::RunControl::Continue
        }
    })
    .map_err(|e| ExecutionError::Failed(e.to_string()))?;

    // 用户在扫描阶段取消时直接结束，避免把部分扫描结果写成“完成输出”。
    if canceled {
        return Err(ExecutionError::Canceled);
    }

    phase(JobPhase::Writing);
    nlr_cli::write_genome_outputs(
        &job.config.input_fasta,
        &job.output_dir,
        &job.prefix,
        job.flank,
        &mut result,
        job.config.locus_filter.as_ref(),
    )
    .map_err(|e| ExecutionError::Failed(e.to_string()))
}

fn execute_protein(
    job: &ProteinJob,
    mut progress: impl FnMut(u64) -> bool,
    mut phase: impl FnMut(JobPhase),
) -> Result<Vec<PathBuf>, ExecutionError> {
    let started = Instant::now();
    let mut canceled = false;
    let result = nlr_cli::run_protein_tracked(&job.config, |p| {
        if progress(p.scanned_fragments) {
            canceled = true;
            nlr_cli::RunControl::Break
        } else {
            nlr_cli::RunControl::Continue
        }
    })
    .map_err(|e| ExecutionError::Failed(e.to_string()))?;

    if canceled {
        return Err(ExecutionError::Canceled);
    }

    let meta = nlr_cli::ProteinRunMeta {
        input_name: job
            .config
            .input_fasta
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown")
            .to_string(),
        threads: job.config.threads,
        elapsed: started.elapsed(),
        lrr_threshold: job.config.motif_accept_p,
        hmm_files: job
            .config
            .hmm_files
            .iter()
            .map(|p| p.display().to_string())
            .collect(),
    };
    phase(JobPhase::Writing);
    nlr_cli::write_protein_outputs(&job.output_dir, &job.prefix, &result, &meta)
        .map_err(|e| ExecutionError::Failed(e.to_string()))
}

fn finalize(state: Arc<JobState>, outcome: Result<Vec<PathBuf>, ExecutionError>) {
    match outcome {
        Ok(paths) => {
            *state.outputs.lock().unwrap() = Some(paths);
            *state.status.lock().unwrap() = JobStatus::Completed;
            *state.phase.lock().unwrap() = JobPhase::Completed;
        }
        Err(ExecutionError::Canceled) => {
            *state.outputs.lock().unwrap() = None;
            *state.status.lock().unwrap() = JobStatus::Canceled;
            *state.phase.lock().unwrap() = JobPhase::Canceled;
        }
        Err(ExecutionError::Failed(message)) => {
            *state.error.lock().unwrap() = Some(message);
            *state.status.lock().unwrap() = JobStatus::Failed;
            *state.phase.lock().unwrap() = JobPhase::Failed;
        }
    }
}

/// 粗估片段总数（与 CLI 的进度条估算口径一致，仅用于前端进度展示）。
fn estimate_total_fragments(path: &Path, fragment_length: usize, overlap: usize) -> u64 {
    let step = fragment_length.saturating_sub(overlap).max(1) as u64;
    let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    let estimated = if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        size.saturating_mul(4)
    } else {
        size
    };
    estimated.saturating_div(step).max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn estimate_total_fragments_is_at_least_one() {
        assert!(estimate_total_fragments(Path::new("does-not-exist.fa"), 20_000, 2_000) >= 1);
    }

    #[test]
    fn service_cancel_returns_false_for_unknown_job() {
        let service = Service::new();
        assert!(!service.cancel(999));
    }

    #[test]
    fn canceled_job_never_exposes_completed_outputs() {
        let state = Arc::new(JobState {
            kind: JobKind::Genome,
            status: Mutex::new(JobStatus::Running),
            phase: Mutex::new(JobPhase::Scanning),
            scanned: AtomicU64::new(10),
            estimated_total: AtomicU64::new(20),
            error: Mutex::new(None),
            outputs: Mutex::new(None),
            cancel: AtomicBool::new(true),
        });

        finalize(state.clone(), Err(ExecutionError::Canceled));

        assert_eq!(*state.status.lock().unwrap(), JobStatus::Canceled);
        assert_eq!(*state.phase.lock().unwrap(), JobPhase::Canceled);
        assert!(state.outputs.lock().unwrap().is_none());
    }
}
