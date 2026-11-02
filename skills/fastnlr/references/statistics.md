# FastNLR 统计口径

- **位点（locus）**：由播种、合并、延伸与冗余吸收得到的 motif 集合，按 `{seq}_nlrN` 命名。
- **完整 NLR**：DNA 模式要求同时含 P-loop（motif_1）与 LRR 类 motif；蛋白模式要求 NB-ARC 与 LRR。
- **结构域类别 `class`**：折叠相邻同类后得到，如 `TIR-NBARC-LRR`、`CC-NBARC-LRR`。
- **家族 `family`**：TNL / RNL / CNL / NB-only / TIR-only / Other；RNL 与 CNL 共享
  `CC-NBARC-LRR` 架构时以 N 端 CC 来源区分。
- **motif 计数**：按 motif id 统计所有通过显著性阈值的命中。
- **每染色体统计**：`motifs`、`nlrs`、`complete_nlrs` 三项。

报告召回/精确率等指标时，需同时说明金标准口径；不同金标准的绝对召回率不可直接比较。
