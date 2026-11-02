# FastNLR MCP 工具

`fastnlr-mcp` 是 stdio MCP 服务，暴露三个工具。

## fastnlr_annotate_genome

输入：`input_path`（必填）、`output_dir`（必填）、`prefix`、`threads`、`motif_accept_p`、
`motif_prelim_p`、`flank`、`flexible_seed`、`require_ploop`、`checkpoint_dir`。
返回：输出文件列表与 `summary.json` 内容。

## fastnlr_annotate_proteins

输入：`input_path`（必填）、`output_dir`（必填）、`prefix`、`threads`、`batch`。
返回：输出文件列表与 `nlr.tsv` 内容。

## fastnlr_plot

输入：`input_dir`（含 `summary.json`，必填）、`output_dir`（必填）、`genome_fai`、
`centromere_bed`。需要本机有 Rscript 及 ggplot2/jsonlite/patchwork。
返回：生成的 PDF 与 HTML 报告路径。
