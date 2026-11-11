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

GFF3 模式（`-i genome.fasta --gff annotation.gff3`）写入
`<prefix>.gff.cds.fasta`、`.gff.proteins.fasta`、`.gff.domains.tsv`、
`.gff.genes.tsv`、`.gff.nlr.tsv`、`.gff.nlr.gff3`、`.gff.skipped.tsv`、
`.gff.summary.txt` 和 `.gff.summary.json`。结构域同时包含蛋白、CDS 和基因组坐标；
坐标均为 1-based、闭区间。

GFF3 模式使用 HMM-only 分类。`domains.tsv` 的 `source` 固定为 `hmm`。
`genes.tsv` 和 JSON 结果中的关键字段：

- `has_lrr_any`：任意位置存在通过 GA 的 LRR；
- `has_cterm_lrr`：至少一个 LRR 位于 NB-ARC 结构域块之后；
- `complete`：存在 NB-ARC 且存在 C 端 LRR；
- `class`：按真实坐标生成的结构域串，NB-ARC 前的 LRR 不进入 `class`；
- `family`：TNL、CNL、RNL、NB-only 或 TIR-only，不代表完整性。

`summary.json` 顶层字段：`sequence_count`、`motif_total`、`nlr_total`、`nlr_complete`、
`per_chromosome`、`motif_counts`、`nlr_type_counts`、`nlr_family_counts`、`loci`。
`loci` 每项含 `name`、`seq`、`start`、`end`、`strand`、`class`、`family`、`complete`、
`has_stop`、`motif_count`；坐标为零基闭区间。
