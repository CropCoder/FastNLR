#!/usr/bin/env Rscript
#
# FastNLR 发表级统计图件（PDF）
# ------------------------------------------------------------
# 功能：
#   读取 FastNLR 输出的 summary.json，生成：
#   - NLR-Stats-Counts.pdf           全局数量与完整度
#   - NLR-Stats-Composition.pdf      类型与家族组成
#   - NLR-Stats-PerChromosome.pdf    各染色体数量统计
#   - NLR-Stats-Lengths.pdf          位点长度分布
#   - NLR-Stats-Motifs.pdf           motif 使用谱
#
# 用法：
#   Rscript R2-NLR-Stats-Report-V1.R <input_dir> <out_dir>
#   input_dir  FastNLR 输出目录（需含 <prefix>.summary.json）
#   out_dir    图件输出目录（不存在会自动创建）
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
  stop("用法: Rscript R2-NLR-Stats-Report-V1.R <input_dir> <out_dir>")
}
input_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- Sys.getenv("FASTNLR_FONT", "Liberation Sans")
pal <- c("#5cc5b2", "#eccc68", "#9c88ff", "#ff6b6b")
theme_pub <- theme_bw(base_size = 9) +
  theme(
    text = element_text(family = font_family),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black")
  )

summary_file <- list.files(input_dir, pattern = "\\.summary\\.json$", full.names = TRUE)
if (length(summary_file) != 1) {
  stop("输入目录中需恰好一个 summary.json：", input_dir)
}
report <- fromJSON(summary_file[1], simplifyVector = FALSE)

# 1. 全局数量与完整度
counts <- data.frame(
  metric = c("Motifs", "NLR loci", "Complete NLR"),
  value = c(report$motif_total, report$nlr_total, report$nlr_complete),
  stringsAsFactors = FALSE
)
p_counts <- ggplot(counts, aes(x = metric, y = value, fill = metric)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = value), vjust = -0.3, size = 3, family = font_family) +
  scale_fill_manual(values = c("Motifs" = pal[1], "NLR loci" = pal[2], "Complete NLR" = pal[3]), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Count", title = "Global counts") +
  theme_pub
ggsave(file.path(out_dir, "NLR-Stats-Counts.pdf"), p_counts,
  width = 90, height = 80, units = "mm", device = cairo_pdf
)

# 2. 类型与家族组成
type_df <- data.frame(
  class = names(report$nlr_type_counts),
  count = unlist(report$nlr_type_counts),
  stringsAsFactors = FALSE
)
type_df <- type_df[order(-type_df$count), ]
type_df$class <- factor(type_df$class, levels = rev(type_df$class))
p_type <- ggplot(type_df, aes(x = class, y = count, fill = class)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = count), hjust = -0.2, size = 3, family = font_family) +
  coord_flip() +
  scale_fill_manual(values = rep(pal, length.out = nrow(type_df)), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "NLR loci", title = "NLR types") +
  theme_pub

family_df <- data.frame(
  family = names(report$nlr_family_counts),
  count = unlist(report$nlr_family_counts),
  stringsAsFactors = FALSE
)
family_df <- family_df[order(-family_df$count), ]
family_df$family <- factor(family_df$family, levels = rev(family_df$family))
p_family <- ggplot(family_df, aes(x = family, y = count, fill = family)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = count), hjust = -0.2, size = 3, family = font_family) +
  coord_flip() +
  scale_fill_manual(values = rep(pal, length.out = nrow(family_df)), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "NLR loci", title = "NLR families") +
  theme_pub

p_comp <- p_type | p_family
ggsave(file.path(out_dir, "NLR-Stats-Composition.pdf"), p_comp,
  width = 180, height = 90, units = "mm", device = cairo_pdf
)

# 3. 各染色体数量统计
chr_df <- do.call(rbind, lapply(names(report$per_chromosome), function(chr) {
  data.frame(
    chr = chr,
    nlrs = report$per_chromosome[[chr]]$nlrs,
    complete = report$per_chromosome[[chr]]$complete_nlrs,
    stringsAsFactors = FALSE
  )
}))
chr_df <- chr_df[order(-chr_df$nlrs), ]
chr_df$chr <- factor(chr_df$chr, levels = chr_df$chr)
chr_long <- reshape(
  chr_df, direction = "long",
  varying = c("nlrs", "complete"),
  v.names = "count", idvar = "chr", timevar = "kind", times = c("All", "Complete")
)
p_chr <- ggplot(chr_long, aes(x = chr, y = count, fill = kind)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.6, family = font_family) +
  scale_fill_manual(values = c("All" = pal[1], "Complete" = pal[3]), name = "Locus") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "NLR loci", title = "NLR loci per chromosome") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out_dir, "NLR-Stats-PerChromosome.pdf"), p_chr,
  width = 180, height = 90, units = "mm", device = cairo_pdf
)

# 4. 位点长度分布
loci <- report$loci
len_df <- data.frame(
  length = sapply(loci, function(x) x$end - x$start + 1),
  family = sapply(loci, function(x) x$family),
  stringsAsFactors = FALSE
)
p_len <- ggplot(len_df, aes(x = length, fill = family)) +
  geom_histogram(bins = 30, colour = "black", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "TNL" = "#5cc5b2", "RNL" = "#9c88ff", "CNL" = "#eccc68",
      "NB-only" = "#b8c2cc", "TIR-only" = "#ff6b6b", "Other" = "#9aa4b0"
    ),
    name = "NLR family", drop = FALSE
  ) +
  scale_x_continuous(labels = function(x) prettyNum(x, big.mark = ",")) +
  labs(x = "Locus length (bp)", y = "Count", title = "NLR locus length distribution") +
  theme_pub
ggsave(file.path(out_dir, "NLR-Stats-Lengths.pdf"), p_len,
  width = 180, height = 90, units = "mm", device = cairo_pdf
)

# 5. motif 使用谱
motif_df <- data.frame(
  motif = sprintf("motif_%d", as.integer(names(report$motif_counts))),
  count = unlist(report$motif_counts),
  stringsAsFactors = FALSE
)
motif_df <- motif_df[order(as.integer(sub("motif_", "", motif_df$motif))), ]
motif_df$motif <- factor(motif_df$motif, levels = motif_df$motif)
p_motif <- ggplot(motif_df, aes(x = motif, y = count)) +
  geom_col(fill = pal[4], colour = "black", linewidth = 0.25) +
  geom_text(aes(label = count), vjust = -0.3, size = 2.4, family = font_family) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Hits", title = "Motif usage") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
ggsave(file.path(out_dir, "NLR-Stats-Motifs.pdf"), p_motif,
  width = 180, height = 90, units = "mm", device = cairo_pdf
)

cat("written statistics PDFs under:", out_dir, "\n")
