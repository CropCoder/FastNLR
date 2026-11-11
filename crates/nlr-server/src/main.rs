//! fastnlr-server — FastNLR 自部署 Web 服务。
//!
//! 以单二进制提供 REST API 与内嵌前端页面，供浏览器与 TBtools 插件使用。
//! 默认只绑定本机回环地址；作业通过 nlr-service 在后台线程执行。

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use axum::extract::{Path as AxumPath, State};
use axum::http::{header, HeaderValue, StatusCode};
use axum::response::{Html, IntoResponse, Response};
use axum::routing::{get, post};
use axum::{Json, Router};
use clap::Parser;
use nlr_service::Service;
use serde::Deserialize;

/// Web 服务命令行参数。
#[derive(Debug, Parser)]
#[command(name = "fastnlr-server", version, about = "FastNLR web server")]
struct Args {
    /// 监听地址（默认仅本机）。
    #[arg(long, default_value = "127.0.0.1")]
    host: String,

    /// 监听端口。
    #[arg(long, default_value_t = 8787)]
    port: u16,

    /// 作业输出根目录（未指定作业输出目录时使用）。
    #[arg(long, default_value = "fastnlr-jobs")]
    output_dir: PathBuf,
}

#[derive(Clone)]
struct AppState {
    service: Service,
    base_output_dir: Arc<PathBuf>,
}

#[derive(Deserialize)]
struct GenomeJobRequest {
    input_path: String,
    #[serde(default)]
    output_dir: Option<String>,
    #[serde(default)]
    prefix: Option<String>,
    #[serde(default)]
    threads: Option<usize>,
    #[serde(default)]
    motif_accept_p: Option<f64>,
    #[serde(default)]
    motif_prelim_p: Option<f64>,
    #[serde(default)]
    flank: Option<u64>,
    #[serde(default)]
    flexible_seed: Option<usize>,
    #[serde(default)]
    no_flexible_seed: bool,
    #[serde(default)]
    require_ploop: bool,
    #[serde(default)]
    mot_file: Option<String>,
    #[serde(default)]
    store_file: Option<String>,
    #[serde(default)]
    checkpoint_dir: Option<String>,
}

#[derive(Deserialize)]
struct AnnotationJobRequest {
    input_path: String,
    gff_path: String,
    #[serde(default)]
    output_dir: Option<String>,
    #[serde(default)]
    prefix: Option<String>,
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

#[derive(serde::Serialize)]
struct JobCreated {
    job_id: u64,
    output_dir: String,
}

#[derive(serde::Serialize)]
struct OutputsResponse {
    outputs: Option<Vec<String>>,
}

#[tokio::main]
async fn main() {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let args = Args::parse();
    let state = AppState {
        service: Service::new(),
        base_output_dir: Arc::new(args.output_dir.clone()),
    };

    let app = Router::new()
        .route("/", get(index))
        .route("/api/health", get(health))
        .route("/api/jobs/genome", post(create_genome_job))
        .route("/api/jobs/annotation", post(create_annotation_job))
        .route("/api/jobs/{id}", get(job_status))
        .route("/api/jobs/{id}/cancel", post(cancel_job))
        .route("/api/jobs/{id}/outputs", get(job_outputs))
        .route("/api/jobs/{id}/files/{filename}", get(download_file))
        .with_state(state);

    let listener = tokio::net::TcpListener::bind((args.host.as_str(), args.port))
        .await
        .unwrap_or_else(|e| panic!("cannot bind {}:{}: {e}", args.host, args.port));
    tracing::info!("fastnlr-server listening on {}:{}", args.host, args.port);
    axum::serve(listener, app).await.unwrap();
}

async fn index() -> Html<&'static str> {
    Html(include_str!("../web/index.html"))
}

async fn health() -> Json<serde_json::Value> {
    Json(serde_json::json!({ "status": "ok" }))
}

async fn create_genome_job(
    State(state): State<AppState>,
    Json(req): Json<GenomeJobRequest>,
) -> Result<Json<JobCreated>, ApiError> {
    let mut config = nlr_service::RunConfig::new(
        PathBuf::from(&req.input_path),
        req.mot_file.map(PathBuf::from),
        req.store_file.map(PathBuf::from),
    );
    if let Some(t) = req.threads {
        config.threads = t;
    }
    if let Some(p) = req.motif_accept_p {
        config.motif_accept_p = p;
    }
    if let Some(p) = req.motif_prelim_p {
        config.motif_prelim_p = p;
    }
    if let Some(cp) = req.checkpoint_dir {
        config.checkpoint_dir = Some(PathBuf::from(cp));
    }
    config.flexible_seed =
        nlr_cli::resolve_flexible_seed(req.flexible_seed, req.no_flexible_seed, req.require_ploop);
    config.validate().map_err(ApiError::bad_request)?;

    let output_dir = resolve_output_dir(&state, req.output_dir.as_deref());
    let prefix = req
        .prefix
        .clone()
        .unwrap_or_else(|| derive_prefix(&req.input_path));
    let job = nlr_service::GenomeJob {
        config,
        output_dir: output_dir.clone(),
        prefix,
        flank: req.flank.unwrap_or(2000),
    };
    let job_id = state.service.submit_genome(job);
    Ok(Json(JobCreated {
        job_id,
        output_dir: output_dir.display().to_string(),
    }))
}

async fn create_annotation_job(
    State(state): State<AppState>,
    Json(req): Json<AnnotationJobRequest>,
) -> Result<Json<JobCreated>, ApiError> {
    if req.motif_accept_p.is_some()
        || req.motif_prelim_p.is_some()
        || req.mot_file.is_some()
        || req.store_file.is_some()
    {
        return Err(ApiError::bad_request(
            "GFF3 annotation mode uses HMM-only classification; PWM fields \
             motif_accept_p, motif_prelim_p, mot_file, and store_file are not supported"
                .to_string(),
        ));
    }
    let mut config = nlr_cli::AnnotationConfig::new(
        PathBuf::from(&req.input_path),
        PathBuf::from(&req.gff_path),
    );
    if let Some(t) = req.threads {
        config.threads = t;
    }
    if let Some(b) = req.batch {
        config.batch = b;
    }
    config.hmm_files = req.hmm_files.iter().map(PathBuf::from).collect();
    config.validate().map_err(ApiError::bad_request)?;

    let output_dir = resolve_output_dir(&state, req.output_dir.as_deref());
    let prefix = req
        .prefix
        .clone()
        .unwrap_or_else(|| derive_prefix(&req.input_path));
    let job = nlr_service::AnnotationJob {
        config,
        output_dir: output_dir.clone(),
        prefix,
    };
    let job_id = state.service.submit_annotation(job);
    Ok(Json(JobCreated {
        job_id,
        output_dir: output_dir.display().to_string(),
    }))
}

async fn job_status(
    State(state): State<AppState>,
    AxumPath(id): AxumPath<u64>,
) -> Result<Json<nlr_service::JobSnapshot>, StatusCode> {
    state
        .service
        .snapshot(id)
        .map(Json)
        .ok_or(StatusCode::NOT_FOUND)
}

async fn cancel_job(
    State(state): State<AppState>,
    AxumPath(id): AxumPath<u64>,
) -> Result<Json<serde_json::Value>, StatusCode> {
    if state.service.cancel(id) {
        Ok(Json(serde_json::json!({ "canceled": true })))
    } else {
        Err(StatusCode::NOT_FOUND)
    }
}

async fn job_outputs(
    State(state): State<AppState>,
    AxumPath(id): AxumPath<u64>,
) -> Result<Json<OutputsResponse>, StatusCode> {
    state
        .service
        .outputs(id)
        .map(|paths| {
            Json(OutputsResponse {
                outputs: Some(paths.iter().map(|p| p.display().to_string()).collect()),
            })
        })
        .map_or_else(|| Err(StatusCode::NOT_FOUND), Ok)
}

async fn download_file(
    State(state): State<AppState>,
    AxumPath((id, filename)): AxumPath<(u64, String)>,
) -> Result<Response, StatusCode> {
    if filename.contains('/') || filename.contains('\\') || filename.contains("..") {
        return Err(StatusCode::BAD_REQUEST);
    }
    let outputs = state.service.outputs(id).ok_or(StatusCode::NOT_FOUND)?;
    let path = outputs
        .iter()
        .find(|p| p.file_name().and_then(|n| n.to_str()) == Some(filename.as_str()))
        .ok_or(StatusCode::NOT_FOUND)?;
    let bytes = tokio::fs::read(path)
        .await
        .map_err(|_| StatusCode::NOT_FOUND)?;
    let mut response = bytes.into_response();
    response.headers_mut().insert(
        header::CONTENT_TYPE,
        HeaderValue::from_static("application/octet-stream"),
    );
    response.headers_mut().insert(
        header::CONTENT_DISPOSITION,
        HeaderValue::from_str(&format!("attachment; filename=\"{filename}\"")).unwrap(),
    );
    Ok(response)
}

struct ApiError {
    status: StatusCode,
    message: String,
}

impl ApiError {
    fn bad_request(message: impl ToString) -> Self {
        ApiError {
            status: StatusCode::BAD_REQUEST,
            message: message.to_string(),
        }
    }
}

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        (
            self.status,
            Json(serde_json::json!({ "error": self.message })),
        )
            .into_response()
    }
}

fn resolve_output_dir(state: &AppState, requested: Option<&str>) -> PathBuf {
    if let Some(dir) = requested {
        let path = PathBuf::from(dir);
        if let Err(e) = std::fs::create_dir_all(&path) {
            tracing::warn!("cannot create output dir {}: {e}", path.display());
        }
        return path;
    }
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_millis();
    let path = state.base_output_dir.join(format!("job-{nonce}"));
    if let Err(e) = std::fs::create_dir_all(&path) {
        tracing::warn!("cannot create output dir {}: {e}", path.display());
    }
    path
}

fn derive_prefix(input: &str) -> String {
    let name = Path::new(input)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("genome");
    let name = name.strip_suffix(".gz").unwrap_or(name);
    name.strip_suffix(".fasta")
        .or_else(|| name.strip_suffix(".fna"))
        .or_else(|| name.strip_suffix(".fas"))
        .or_else(|| name.strip_suffix(".fa"))
        .unwrap_or(name)
        .to_string()
}
