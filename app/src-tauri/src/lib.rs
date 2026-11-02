//! FastNLR 桌面版后端：通过 nlr-service 暴露基因组/蛋白注释命令。

use nlr_service::{GenomeJob, ProteinJob};
use serde::Serialize;

#[derive(Serialize)]
pub struct RunOutputs {
    pub outputs: Vec<String>,
}

#[derive(serde::Deserialize)]
pub struct GenomeArgs {
    pub input_path: String,
    pub output_dir: String,
    pub prefix: Option<String>,
    pub threads: Option<usize>,
    pub flank: Option<u64>,
}

#[derive(serde::Deserialize)]
pub struct ProteinArgs {
    pub input_path: String,
    pub output_dir: String,
    pub prefix: Option<String>,
    pub threads: Option<usize>,
}

#[tauri::command]
async fn run_genome(args: GenomeArgs) -> Result<RunOutputs, String> {
    tauri::async_runtime::spawn_blocking(move || {
        let mut config = nlr_service::RunConfig::new(
            std::path::PathBuf::from(&args.input_path),
            None,
            None,
        );
        if let Some(t) = args.threads {
            config.threads = t;
        }
        config.flexible_seed = nlr_cli::resolve_flexible_seed(None, false, false);
        let prefix = args
            .prefix
            .unwrap_or_else(|| derive_prefix(&args.input_path));
        let job = GenomeJob {
            config,
            output_dir: std::path::PathBuf::from(&args.output_dir),
            prefix,
            flank: args.flank.unwrap_or(2000),
        };
        nlr_service::run_genome_blocking(&job).map(|outputs| RunOutputs {
            outputs: outputs.iter().map(|p| p.display().to_string()).collect(),
        })
    })
    .await
    .map_err(|e| e.to_string())?
}

#[tauri::command]
async fn run_protein(args: ProteinArgs) -> Result<RunOutputs, String> {
    tauri::async_runtime::spawn_blocking(move || {
        let mut config = nlr_cli::ProteinConfig::new(std::path::PathBuf::from(&args.input_path));
        if let Some(t) = args.threads {
            config.threads = t;
        }
        let prefix = args
            .prefix
            .unwrap_or_else(|| derive_prefix(&args.input_path));
        let job = ProteinJob {
            config,
            output_dir: std::path::PathBuf::from(&args.output_dir),
            prefix,
        };
        nlr_service::run_protein_blocking(&job).map(|outputs| RunOutputs {
            outputs: outputs.iter().map(|p| p.display().to_string()).collect(),
        })
    })
    .await
    .map_err(|e| e.to_string())?
}

fn derive_prefix(input: &str) -> String {
    let name = std::path::Path::new(input)
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

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .invoke_handler(tauri::generate_handler![run_genome, run_protein])
        .run(tauri::generate_context!())
        .expect("error while running FastNLR desktop app");
}
