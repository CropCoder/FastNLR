//! fastnlr-mcp — 面向 Codex/Claude Code 等智能体的 MCP stdio 服务。
//!
//! 以新行分隔的 JSON-RPC 2.0 运行在标准输入/输出上，暴露 FastNLR 的基因组注释、
//! 蛋白注释与 R 可视化工具。第一阶段工具为同步阻塞调用，长任务通过检查点续跑配合。

use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::process::Command;

use serde_json::{json, Value};

const PROTOCOL_VERSION: &str = "2025-03-26";
const TOOL_VERSION: &str = env!("CARGO_PKG_VERSION");

/// 内嵌的 R 可视化脚本（避免要求调用方额外下载脚本文件）。
const R0: &str = include_str!("../../../visualization/R0-Run-FastNLR-Visualization-V1.R");
const R1: &str = include_str!("../../../visualization/R1-NLR-Chromosome-Map-V1.R");
const R2: &str = include_str!("../../../visualization/R2-NLR-Stats-Report-V1.R");

fn main() {
    let stdin = std::io::stdin();
    let mut stdout = std::io::stdout();
    for line in stdin.lock().lines() {
        let Ok(line) = line else { break };
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let Ok(request) = serde_json::from_str::<Value>(line) else {
            continue;
        };
        let Some(method) = request.get("method").and_then(Value::as_str) else {
            continue;
        };
        let id = request.get("id").cloned();
        let params = request.get("params").cloned().unwrap_or(Value::Null);

        let response = match dispatch(method, params) {
            Ok(result) => json!({ "jsonrpc": "2.0", "id": id, "result": result }),
            Err(message) => json!({
                "jsonrpc": "2.0",
                "id": id,
                "error": { "code": -32000, "message": message }
            }),
        };
        // 通知（无 id）不回复。
        if id.is_some() {
            writeln!(stdout, "{response}").ok();
            stdout.flush().ok();
        }
    }
}

fn dispatch(method: &str, params: Value) -> Result<Value, String> {
    match method {
        "initialize" => Ok(json!({
            "protocolVersion": PROTOCOL_VERSION,
            "capabilities": { "tools": { "listChanged": false } },
            "serverInfo": { "name": "fastnlr-mcp", "version": TOOL_VERSION }
        })),
        "ping" => Ok(json!({})),
        "tools/list" => Ok(tools_list()),
        "tools/call" => call_tool(params),
        _ => Err(format!("unknown method: {method}")),
    }
}

fn tools_list() -> Value {
    json!({
        "tools": [
            {
                "name": "fastnlr_annotate_genome",
                "description": "Annotate NLR loci in a plant genome FASTA and return run statistics plus output files.",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "input_path": { "type": "string", "description": "Genome FASTA path (gzip allowed)" },
                        "output_dir": { "type": "string" },
                        "prefix": { "type": "string" },
                        "threads": { "type": "integer", "minimum": 1 },
                        "motif_accept_p": { "type": "number" },
                        "motif_prelim_p": { "type": "number" },
                        "flank": { "type": "integer" },
                        "flexible_seed": { "type": "integer" },
                        "require_ploop": { "type": "boolean" },
                        "checkpoint_dir": { "type": "string" }
                    },
                    "required": ["input_path", "output_dir"]
                }
            },
            {
                "name": "fastnlr_annotate_proteins",
                "description": "Classify NLR domain architectures in a protein FASTA with bundled HMMER models and PWM LRR.",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "input_path": { "type": "string" },
                        "output_dir": { "type": "string" },
                        "prefix": { "type": "string" },
                        "threads": { "type": "integer", "minimum": 1 },
                        "batch": { "type": "integer", "minimum": 1 }
                    },
                    "required": ["input_path", "output_dir"]
                }
            },
            {
                "name": "fastnlr_plot",
                "description": "Generate publication-style PDF figures and an HTML report from a completed FastNLR run (requires Rscript).",
                "inputSchema": {
                    "type": "object",
                    "properties": {
                        "input_dir": { "type": "string", "description": "FastNLR output directory containing summary.json" },
                        "output_dir": { "type": "string" },
                        "genome_fai": { "type": "string" },
                        "centromere_bed": { "type": "string" }
                    },
                    "required": ["input_dir", "output_dir"]
                }
            }
        ]
    })
}

fn call_tool(params: Value) -> Result<Value, String> {
    let name = params
        .get("name")
        .and_then(Value::as_str)
        .ok_or("tools/call missing name")?;
    let arguments = params.get("arguments").cloned().unwrap_or(json!({}));
    let text = match name {
        "fastnlr_annotate_genome" => annotate_genome(&arguments)?,
        "fastnlr_annotate_proteins" => annotate_proteins(&arguments)?,
        "fastnlr_plot" => plot(&arguments)?,
        _ => return Err(format!("unknown tool: {name}")),
    };
    Ok(json!({
        "content": [{ "type": "text", "text": text }],
        "isError": false
    }))
}

fn annotate_genome(arguments: &Value) -> Result<String, String> {
    let input = required_str(arguments, "input_path")?;
    let output_dir = required_str(arguments, "output_dir")?;
    let mut config = nlr_service::RunConfig::new(
        PathBuf::from(input),
        optional_path(arguments, "mot_file"),
        optional_path(arguments, "store_file"),
    );
    if let Some(t) = get_usize(arguments, "threads") {
        config.threads = t;
    }
    if let Some(p) = get_f64(arguments, "motif_accept_p") {
        config.motif_accept_p = p;
    }
    if let Some(p) = get_f64(arguments, "motif_prelim_p") {
        config.motif_prelim_p = p;
    }
    if let Some(cp) = optional_path(arguments, "checkpoint_dir") {
        config.checkpoint_dir = Some(cp);
    }
    config.flexible_seed = nlr_cli::resolve_flexible_seed(
        get_usize(arguments, "flexible_seed"),
        false,
        get_bool(arguments, "require_ploop"),
    );
    config.validate().map_err(|e| e.to_string())?;

    let job = nlr_service::GenomeJob {
        config,
        output_dir: PathBuf::from(output_dir),
        prefix: get_str(arguments, "prefix")
            .map(str::to_string)
            .unwrap_or_else(|| derive_prefix(input)),
        flank: get_u64(arguments, "flank").unwrap_or(2000),
    };
    let outputs = nlr_service::run_genome_blocking(&job)?;
    read_summary_json(&job.output_dir).map(|summary| {
        format!(
            "FastNLR genome annotation completed.\noutputs:\n{}\nsummary:\n{}",
            outputs
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join("\n"),
            summary
        )
    })
}

fn annotate_proteins(arguments: &Value) -> Result<String, String> {
    let input = required_str(arguments, "input_path")?;
    let output_dir = required_str(arguments, "output_dir")?;
    let mut config = nlr_cli::ProteinConfig::new(PathBuf::from(input));
    if let Some(t) = get_usize(arguments, "threads") {
        config.threads = t;
    }
    if let Some(b) = get_usize(arguments, "batch") {
        config.batch = b;
    }
    config.validate().map_err(|e| e.to_string())?;

    let job = nlr_service::ProteinJob {
        config,
        output_dir: PathBuf::from(output_dir),
        prefix: get_str(arguments, "prefix")
            .map(str::to_string)
            .unwrap_or_else(|| derive_prefix(input)),
    };
    let outputs = nlr_service::run_protein_blocking(&job)?;
    let nlr_tsv = job.output_dir.join(format!("{}.nlr.tsv", job.prefix));
    let summary = std::fs::read_to_string(nlr_tsv).unwrap_or_else(|_| String::new());
    Ok(format!(
        "FastNLR protein annotation completed.\noutputs:\n{}\nNLR candidates (TSV):\n{}",
        outputs
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join("\n"),
        summary
    ))
}

fn plot(arguments: &Value) -> Result<String, String> {
    let input_dir = required_str(arguments, "input_dir")?;
    let output_dir = required_str(arguments, "output_dir")?;
    let tmp = std::env::temp_dir().join(format!("fastnlr-mcp-r-{}", std::process::id()));
    std::fs::create_dir_all(&tmp).map_err(|e| e.to_string())?;
    std::fs::write(tmp.join("R0-Run-FastNLR-Visualization-V1.R"), R0).map_err(|e| e.to_string())?;
    std::fs::write(tmp.join("R1-NLR-Chromosome-Map-V1.R"), R1).map_err(|e| e.to_string())?;
    std::fs::write(tmp.join("R2-NLR-Stats-Report-V1.R"), R2).map_err(|e| e.to_string())?;

    let mut args: Vec<String> = vec![
        tmp.join("R0-Run-FastNLR-Visualization-V1.R")
            .display()
            .to_string(),
        input_dir.to_string(),
        output_dir.to_string(),
    ];
    if let Some(fai) = optional_path(arguments, "genome_fai") {
        args.push(fai.display().to_string());
        if let Some(bed) = optional_path(arguments, "centromere_bed") {
            args.push(bed.display().to_string());
        }
    }

    let output = Command::new("Rscript")
        .args(&args)
        .output()
        .map_err(|e| format!("Rscript not available: {e}"))?;
    if !output.status.success() {
        return Err(format!(
            "plot failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    Ok(format!(
        "plots written to {output_dir}:\n{}",
        String::from_utf8_lossy(&output.stdout)
    ))
}

fn read_summary_json(output_dir: &std::path::Path) -> Result<String, String> {
    let entries = std::fs::read_dir(output_dir).map_err(|e| e.to_string())?;
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().and_then(|e| e.to_str()) == Some("json") {
            return std::fs::read_to_string(&path).map_err(|e| e.to_string());
        }
    }
    Err("summary.json not found in output directory".to_string())
}

fn required_str<'a>(arguments: &'a Value, key: &str) -> Result<&'a str, String> {
    arguments
        .get(key)
        .and_then(Value::as_str)
        .filter(|s| !s.trim().is_empty())
        .ok_or_else(|| format!("missing required argument: {key}"))
}

fn get_str<'a>(arguments: &'a Value, key: &str) -> Option<&'a str> {
    arguments.get(key).and_then(Value::as_str)
}

fn optional_path(arguments: &Value, key: &str) -> Option<PathBuf> {
    get_str(arguments, key)
        .filter(|s| !s.trim().is_empty())
        .map(PathBuf::from)
}

fn get_u64(arguments: &Value, key: &str) -> Option<u64> {
    arguments.get(key).and_then(Value::as_u64)
}

fn get_usize(arguments: &Value, key: &str) -> Option<usize> {
    arguments
        .get(key)
        .and_then(Value::as_u64)
        .map(|v| v as usize)
}

fn get_f64(arguments: &Value, key: &str) -> Option<f64> {
    arguments.get(key).and_then(Value::as_f64)
}

fn get_bool(arguments: &Value, key: &str) -> bool {
    arguments.get(key).and_then(Value::as_bool).unwrap_or(false)
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
