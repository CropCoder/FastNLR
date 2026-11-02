#!/usr/bin/env Rscript
#
# FastNLR 可视化总入口（投稿级 PDF + HTML 报告）
# ------------------------------------------------------------
# 功能：
#   依次调用 R1（染色体全景图）与 R2（统计图件），并生成 HTML 报告索引。
#
# 用法：
#   Rscript R0-Run-FastNLR-Visualization-V1.R <input_dir> <out_dir> [<genome.fai>] [<centromere.bed>]
#   input_dir      FastNLR 输出目录（需含 <prefix>.summary.json）
#   out_dir        图件与报告输出目录
#   genome.fai     可选，染色体长度索引（samtools faidx 生成）
#   centromere.bed 可选，着丝粒/缺口 BED
#
# Author: Jiwen Zhao (https://github.com/CropCoder)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法: Rscript R0-Run-FastNLR-Visualization-V1.R <input_dir> <out_dir> [<genome.fai>] [<centromere.bed>]")
}
input_dir <- args[1]
out_dir <- args[2]
extra_args <- if (length(args) >= 3) args[3:length(args)] else character(0)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (!nzchar(script_dir)) {
  script_dir <- "."
}

r1 <- file.path(script_dir, "R1-NLR-Chromosome-Map-V1.R")
r2 <- file.path(script_dir, "R2-NLR-Stats-Report-V1.R")

system2("Rscript", c(r1, input_dir, out_dir, extra_args))
system2("Rscript", c(r2, input_dir, out_dir))

summary_file <- list.files(input_dir, pattern = "\\.summary\\.json$", full.names = TRUE)
report <- if (length(summary_file) == 1) jsonlite::fromJSON(summary_file[1], simplifyVector = FALSE) else list()
pdfs <- sort(list.files(out_dir, pattern = "\\.pdf$", full.names = TRUE))

html <- c(
  "<!doctype html><html><head><meta charset='utf-8'><title>FastNLR report</title></head><body>",
  "<h1>FastNLR report</h1>",
  sprintf(
    "<p>NLR loci: %s | Complete: %s | Motifs: %s</p>",
    ifelse(is.null(report$nlr_total), 0, report$nlr_total),
    ifelse(is.null(report$nlr_complete), 0, report$nlr_complete),
    ifelse(is.null(report$motif_total), 0, report$motif_total)
  ),
  "<ul>",
  sprintf("<li><a href='%s'>%s</a></li>", basename(pdfs), basename(pdfs)),
  "</ul></body></html>"
)
writeLines(html, file.path(out_dir, "NLR-Report.html"))
cat("written:", file.path(out_dir, "NLR-Report.html"), "\n")
