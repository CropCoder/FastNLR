# FastNLR 输出格式

基因组模式（输入 `genome.fasta`，前缀 `genome`）默认写入：

| 文件 | 内容 |
|---|---|
| `genome.nlr.txt` | 位点表：seqname、name、class、start、end、strand |
| `genome.nlr.gff` | GFF2 位点记录，feature `NBSLRR` |
| `genome.nlr.bed` | BED12，每 motif 一个 block，颜色区分完整/部分/含终止子 |
| `genome.motifs.bed` | 逐 motif 命中，score 列为 p 值 |
| `genome.nbarc.fasta` | P-loop 锚定的 NB-ARC 多重比对 |
| `genome.loci.fasta` | 位点侧翼序列（±flank），负链反向互补 |
| `genome.stats.tsv` | 全局/染色体/motif/家族统计 |
| `genome.summary.json` | 机器可读汇总 |
| `genome.summary.txt` | 人类可读运行摘要 |

蛋白模式写入 `domains.tsv`、`nlr.tsv`、`summary.txt`，不产生坐标类文件。

`summary.json` 顶层字段：`sequence_count`、`motif_total`、`nlr_total`、`nlr_complete`、
`per_chromosome`、`motif_counts`、`nlr_type_counts`、`nlr_family_counts`、`loci`。
`loci` 每项含 `name`、`seq`、`start`、`end`、`strand`、`class`、`family`、`complete`、
`has_stop`、`motif_count`；坐标为零基闭区间。
