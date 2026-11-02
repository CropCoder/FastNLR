#!/usr/bin/env Rscript
#
# FastNLR 染色体全景图（投稿级 PDF）
# ------------------------------------------------------------
# 功能：
#   读取 FastNLR 输出的 summary.json，绘制染色体全景图：
#   - 每条染色体按长度绘制为水平条（长度来自 .fai，缺省时回退到位点最大坐标）
#   - 可选标出着丝粒/缺口区间（外部 BED，0-based 半开区间）
#   - 按生物学家族（TNL/RNL/CNL/NB-only/TIR-only/Other）着色展示 NLR 位点
#   - 正链画在染色体线上方、负链画在下方
#
# 用法：
#   Rscript R1-NLR-Chromosome-Map-V1.R <input_dir> <out_dir> [<genome.fai>] [<centromere.bed>]
#   input_dir      FastNLR 输出目录（需含 <prefix>.summary.json）
#   out_dir        图件输出目录（不存在会自动创建）
#   genome.fai     可选，samtools faidx 生成的染色体索引，用于精确长度
#   centromere.bed 可选，三列 BED：chrom start end
#
# 说明：
#   - 图中不出现中文或特殊字符，字体使用 Arial 的度量兼容替代 Liberation Sans；
#     若系统装有 Arial，可设置环境变量 FASTNLR_FONT=Arial。
#   - 主题色遵循项目规范：#5cc5b2 #eccc68 #9c88ff #ff6b6b。
#
# Author: Jiwen Zhao (https://github.com/CropCoder)

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("用法: Rscript R1-NLR-Chromosome-Map-V1.R <input_dir> <out_dir> [<genome.fai>] [<centromere.bed>]")
}
input_dir <- args[1]
out_dir <- args[2]
fai_path <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL
centromere_path <- if (length(args) >= 4 && nzchar(args[4])) args[4] else NULL

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- Sys.getenv("FASTNLR_FONT", "Liberation Sans")
pal_family <- c(
  "TNL" = "#5cc5b2",
  "RNL" = "#9c88ff",
  "CNL" = "#eccc68",
  "NB-only" = "#b8c2cc",
  "TIR-only" = "#ff6b6b",
  "Other" = "#9aa4b0"
)

summary_file <- list.files(input_dir, pattern = "\\.summary\\.json$", full.names = TRUE)
if (length(summary_file) != 1) {
  stop("输入目录中需恰好一个 summary.json：", input_dir)
}
report <- fromJSON(summary_file[1], simplifyVector = FALSE)

chromosomes <- names(report$per_chromosome)
loci <- report$loci

# 染色体长度：优先 .fai，否则回退到位点最大坐标。
if (!is.null(fai_path)) {
  fai <- read.table(fai_path, sep = "\t", quote = "", stringsAsFactors = FALSE)
  chr_len <- setNames(as.numeric(fai[[2]]), fai[[1]])
  chr_order <- intersect(names(chr_len), chromosomes)
  if (length(chr_order) == 0) {
    chr_order <- chromosomes
  }
} else {
  chr_len <- NULL
  chr_order <- chromosomes
}

max_end_by_chr <- tapply(loci, sapply(loci, function(x) x$seq), function(grp) {
  max(sapply(grp, function(x) x$end))
})

chr_tab <- do.call(rbind, lapply(chr_order, function(chr) {
  if (!is.null(chr_len) && chr %in% names(chr_len)) {
    len <- chr_len[[chr]]
  } else {
    len <- max_end_by_chr[[chr]] * 1.05
  }
  data.frame(chr = chr, len_bp = len, stringsAsFactors = FALSE)
}))
chr_tab <- chr_tab[order(-chr_tab$len_bp), ]
chr_tab$y <- seq_len(nrow(chr_tab))

map_df <- chr_tab
loci_df <- data.frame(
  seq = sapply(loci, function(x) x$seq),
  mid = sapply(loci, function(x) (x$start + x$end) / 2),
  strand = sapply(loci, function(x) x$strand),
  family = sapply(loci, function(x) x$family),
  stringsAsFactors = FALSE
)
loci_df <- loci_df[loci_df$seq %in% chr_tab$chr, ]
loci_df$y <- chr_tab$y[match(loci_df$seq, chr_tab$chr)]
loci_df$y0 <- ifelse(loci_df$strand == "+", loci_df$y + 0.08, loci_df$y - 0.08)
loci_df$y1 <- ifelse(loci_df$strand == "+", loci_df$y + 0.34, loci_df$y - 0.34)

p <- ggplot() +
  geom_segment(
    data = map_df,
    aes(x = 0, xend = len_bp / 1e6, y = y, yend = y),
    linewidth = 2.6, colour = "grey87", lineend = "round"
  )

if (!is.null(centromere_path)) {
  cent <- read.table(centromere_path, sep = "\t", quote = "", stringsAsFactors = FALSE)
  colnames(cent) <- c("chr", "start", "end")[seq_len(ncol(cent))]
  cent <- cent[cent$chr %in% chr_tab$chr, ]
  if (nrow(cent) > 0) {
    cent$y <- chr_tab$y[match(cent$chr, chr_tab$chr)]
    p <- p + geom_rect(
      data = cent,
      aes(xmin = start / 1e6, xmax = end / 1e6, ymin = y - 0.12, ymax = y + 0.12),
      fill = "grey45", alpha = 0.75
    )
  }
}

p <- p +
  geom_segment(
    data = loci_df,
    aes(x = mid / 1e6, xend = mid / 1e6, y = y0, yend = y1, colour = family),
    linewidth = 0.55, lineend = "round"
  ) +
  scale_colour_manual(values = pal_family, name = "NLR family", drop = FALSE) +
  scale_y_continuous(breaks = chr_tab$y, labels = chr_tab$chr, expand = expansion(mult = c(0.04, 0.04))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    x = "Position (Mb)", y = NULL,
    title = "NLR loci across chromosomes",
    subtitle = "Upper/lower ticks: forward and reverse strand"
  ) +
  theme_bw(base_size = 9) +
  theme(
    text = element_text(family = font_family),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black")
  )

out_pdf <- file.path(out_dir, "NLR-Chromosome-Map.pdf")
ggsave(out_pdf, p, width = 180, height = max(80, 18 * nrow(chr_tab)), units = "mm", device = cairo_pdf)
cat("written:", out_pdf, "\n")
