---
name: fastnlr
description: Install and run FastNLR to annotate NLR immune-receptor loci from a plant genome FASTA or classify NLR domains from a protein FASTA, then interpret the results and generate publication-style figures. Use when the user asks to annotate NLR genes, scan a genome for NLR loci, classify NLR proteins, or produce FastNLR reports/plots.
---

# FastNLR

FastNLR 是一个高速、自包含的植物 NLR 免疫受体位点注释工具。它扫描基因组六框翻译后的
氨基酸 motif，并组装出完整或部分 NLR 位点；也支持对蛋白 FASTA 做结构域分类。

## 安装

优先使用预编译二进制或 Bioconda，避免不必要的源码编译。

```bash
# 方式一：Bioconda（推荐）
conda install bioconda::fastnlr

# 方式二：本项目脚本（Linux/macOS）
bash scripts/install_fastnlr.sh

# 方式三：源码编译（需要 Rust 1.85+）
git clone https://github.com/CropCoder/FastNLR.git
cd FastNLR && cargo build --release
```

Windows 用户可运行 `scripts/install_fastnlr.ps1`，或从 GitHub Releases 下载对应
`fastnlr-<target>.zip`。

## 基因组模式

```bash
fastnlr -i genome.fasta -o results -t 8
```

常用参数：

- `-i <fasta>`：基因组 FASTA（支持 gzip），必填。
- `-o <dir>`：输出目录，必填。
- `-t <n>`：线程数。
- `--checkpoint <dir>`：断点续跑，跳过已完成的扫描。
- `--motif-accept-p` / `--motif-prelim-p`：放宽可提高对分化 NLR 的召回，但会显著增加耗时。
- `--flexible-seed <n>` 默认开启（n=3）；`--no-flexible-seed` 回退到精确种子组合。

输出文件与字段含义见 [references/output-format.md](references/output-format.md)。

## 蛋白模式

```bash
fastnlr --protein -i proteins.fa -o protein_results
```

该模式用内置 Pfam HMM（NB-ARC、TIR、Rx_N、RPW8、R1-type）加 PWM LRR 判定结构域架构，
输出 `<prefix>.domains.tsv`、`<prefix>.nlr.tsv` 和 `<prefix>.summary.txt`。

## 结果解读

- 结构域类别串（如 `CC-NBARC-LRR`）由 `class` 字段给出。
- 生物学家族 `family` 区分 TNL / RNL / CNL / NB-only / TIR-only / Other。
- `complete` 表示含 NB-ARC 与 LRR（DNA 模式另要求 P-loop）。
- 机器可读汇总在 `<prefix>.summary.json`，统计口径见
  [references/statistics.md](references/statistics.md)。

## 发表级可视化

需要 R 与 `ggplot2`、`jsonlite`、`patchwork`。仓库的 `visualization/` 提供：

```bash
Rscript visualization/R0-Run-FastNLR-Visualization-V1.R \
  results results_viz [genome.fai] [centromere.bed]
```

生成染色体全景图、数量/组成/长度/motif 等 PDF 与 HTML 报告。着丝粒/缺口来自可选 BED。

