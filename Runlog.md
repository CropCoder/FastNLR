# FastNLR Run Log

## 2026-09-21

### Performance benchmark

- Built release binary with `cargo build --release --locked`.
- Created benchmark workflow:
  - `benchmark/Script1-Run-Performance-Benchmark-V1.R`
  - `benchmark/Script2-Plot-Performance-Benchmark-V1.R`
- Benchmark designs:
  - all species, 8 threads, default parameters;
  - representative species, 1/4/8/16 threads;
  - representative species, default/sensitive/relaxed parameter modes.
- Metrics: wall time, maximum RSS, CPU percent, NLR count, motif count.

### Performance benchmark results

- Completed 23 out of 33 planned scenarios.
- Completed species default runs for 15 species; large runs for
  `Gossypium_barbadense`, `Hordeum_vulgare`, and `Zea_mays` were not
  finished during this session and were omitted from plotting.
- Thread scaling and parameter sensitivity are available for the
  completed representative subset (mainly `Cucumis_sativus` and
  `Glycine_max`).
- Raw metrics: `benchmark/metrics/benchmark_raw.csv`.
- Figures: `benchmark/figures/Fig1-*.pdf` through `Fig6-*.pdf`.

### Largest-genome thread scaling

- Tested `Zanthoxylum_bungeanum` with 256, 128, 64, 32, 16, 8, and 4 threads.
- Wall time decreased with more threads, while peak RSS remained around 4.5 GB.
- Threads 1 and 2 were not completed in this session; their rows remain `NA`
  in `benchmark/metrics/largest_thread_scaling_raw.csv`.
- Figure: `benchmark/figures/Fig7-largest-genome-thread-scaling.pdf`.

## 2026-09-22

### GoldNLR / RefPlantNLR 阳性数据集准确度评测

- 数据集：RefPlantNLR 已发表阳性集 `Supplemental_dataset_2_RefPlantNLR_v.20200528_415_CDS.fasta`
  （400 条 NLR CDS）+ `Table_S1_RefPlantNLR_v.20200528_415_Description.xlsx`（415 行，其中 400 行有 CDS）。
- 脚本：
  - `benchmark/Script5-Run-GoldNLR-Accuracy-Eval-V1.R`（运行 + 评估）
  - `benchmark/Script6-Plot-GoldNLR-Accuracy-V1.R`（绘图）
- 配置：default（内置 20 motif）与 extended（`data/motif_library_rnl_tir/` 扩展 28 motif）两档。
- 结果：
  - 基因级召回率：default 393/400（98.25%），extended 398/400（99.50%）。
  - RNL（helper NLR）召回：default 3/8，extended 8/8（扩展库补全 RNL 召回）。
  - 两档均漏检 2 个：ADR1-L3（NL）、Wed（CNL）。
  - 结构域 N 端一致性（extended，三类 CC/TIR/NB-only）：91.88%（BED-CC 型 4 个另归 Other，未纳入）。
  - LRR 检出（extended）：灵敏度 99.74%、特异度 100%。
  - 过度切分（>1 locus）：default 3 个基因、extended 28 个基因，且扩展档多 locus 全为 TNL/TNLN/TNOTNL 类——
    扩展库新增 TIR-only 种子导致 TIR 域与 NB-ARC/LRR 被切分为多个位点，属以 TNL 过度切分换取 RNL 召回。
- 产物：
  - `benchmark/metrics/goldnlr_gene_table.csv` / `goldnlr_summary.csv` / `goldnlr_confusion.csv`
  - `benchmark/figures/Fig8~Fig10-goldnlr-*.pdf`

## 2026-09-23

### 柔性播种（category/rank-based seeding）升级方案定稿

- 定位问题：内置 `is_seed` 采用「整串精确相等」，导致 RNL 类命中串 `[6,4,10,3,2,...]`
  无法播种；RefPlantNLR 阳性集上 default 档 RNL 召回仅 3/8。
- 方案：把「精确 ID 串匹配」升级为「NB-ARC 按 rank 单调 + 数量达标」的柔性播种，
  并同步放宽 `has_signature` 预过滤。
- 改动文件：`crates/nlr-core/src/signature_def.rs`、`crates/nlr-assemble/src/lib.rs`、
  `crates/nlr-cli/src/lib.rs`、`crates/nlr-cli/src/main.rs`。
- 新 CLI：`--flexible-seed <n>`（默认 3）与 `--require-ploop`（默认 false，因目标 RNL
  常缺 P-loop motif_1）。
- 完整设计见 `docs/flexible-seeding-upgrade.md`；默认保持 opt-in，验证后再考虑设为默认。

## 2026-09-23（续）

### 柔性播种落地实现（category/rank-based seeding）

- 按 `docs/flexible-seeding-upgrade.md` 落地代码，共改 5 个文件（+117 / -7）：
  - `crates/nlr-core/src/signature_def.rs`：新增 `FlexibleSeedConfig`、`is_seed_flexible()`、
    `has_signature_flexible()`、`with_flexible_seed()`，`AnnotatorSignatureDefinition` 加
    `flexible_seed` 字段。
  - `crates/nlr-assemble/src/lib.rs`：`find_seeds` 播种判定改为 flexible → relaxed → 精确匹配三态。
  - `crates/nlr-cli/src/lib.rs`：`RunConfig` 加 `flexible_seed` 字段并注入 `signature_def`；
    `scan_one_fragment` 预过滤按需改用 `has_signature_flexible`。
  - `crates/nlr-cli/src/main.rs`：新增 `--flexible-seed <n>` 与 `--require-ploop`。
  - `crates/nlr-core/tests/core.rs`：新增 2 个单元测试。
- 测试：`cargo test --workspace --locked` 全部通过（含新增 `flexible_seed_rank_monotonic`、
  `flexible_signature_filter`）。
- 功能验证（RefPlantNLR 400 条 CDS，`--flexible-seed 3`，默认 20 motif）：
  - 基因级召回 398/400（99.5%），RNL 召回 3/8 → 8/8，过度切分仅 3 个基因（扩展库为 28）。
  - 结论：柔性播种以 rank 单调 + 最少 NB-ARC 数播种，补全 RNL 且不引入 TIR-only 过度切分。
- 说明：本地 rustfmt 1.8.0 与仓库既有格式版本不一致，`cargo fmt --check` 全库报既有漂移；
  为避免改动无关代码，未执行 `cargo fmt`。

## 2026-09-23（再续）

### 输出参数重构：去掉 -p、-o 统一为输出目录、默认全量输出

- 目标：`-o <dir>` 直接指定输出目录（不存在自动创建），去掉 `-p`，所有结果写入该目录，
  文件名前缀取输入基因组文件名（去扩展名）；默认输出全部结果（序列/表格/图/摘要）。
- 改动（`crates/nlr-cli/src/main.rs`）：
  - 移除 `-p/--output-prefix`，以及 `-g/-b/-m/-a/-f/-c/--stats/--plot/--summary` 等可选输出参数。
  - 新增 `--flank <bp>`（默认 2000）控制位点序列提取的侧翼长度。
  - 前缀由 `genome_prefix()` 从 `-i` 文件名推导（处理 .fa/.fasta/.fna/.fas/.gz）。
  - `write_outputs` 统一输出：`{prefix}.nlr.txt/.nlr.gff/.nlr.bed/.motifs.bed/.nbarc.fasta/.loci.fasta/.stats.tsv/.summary.txt` + `{prefix}.plots/`（PNG）。
  - `print_run_summary` 改为 `build_summary`，摘要写 `.summary.txt` 并打印 stdout；补充 per-chromosome 表。
  - 顺带移除进度条结束后的 2 秒 sleep。
- 验证：`cargo test --workspace --locked` 全部通过、无 warning；功能验证 `-o /tmp/newcli_out` 正确生成
  以 `test_genome` 为前缀的全部 8 个文件 + plots 目录。
- 遗留：README 与 benchmark/Script5、Script6 仍使用旧 `-p` / `-o out.txt` 用法，需同步适配。

## 2026-09-24

### 修复绘图空白与坐标问题（nlr-plot）

- 根因：plotters 的 ab_glyph 后端默认不内置任何字体，`("sans-serif", N)` 文字绘制报
  `FontUnavailable`，导致两张 PNG 全部空白（纯白，md5 相同）。
- 修复：
  - 内嵌 `crates/nlr-plot/assets/LiberationSans-Regular.ttf`（SIL OFL 1.1，与 GPL-3.0 兼容），
    在 `nlr-plot` 中用 `register_font` + `Once` 注册为 "sans-serif"。
  - motif 计数图：x 轴由硬编码 `1..21` 改为按实际 motif id 动态范围，并改用
    `into_segmented()` 使柱居中对齐（支持扩展库 21-28）。
  - 染色体 NLR 图：x 轴改用 `into_segmented()` + `x_label_formatter` 显示染色体名，
    替换原先的纯数字索引。
- 验证：`cargo test --workspace --locked` 全部通过、无 warning；功能验证两张 PNG 非白像素
  比例分别 0.19 / 0.12，且 md5 不同（两张图内容不同）。

## 2026-09-24（再续）

### 调整运行摘要（summary）输出格式

- 删除 per-chromosome summary 段落。
- NLR types 列顺序改为「数量在前、类型在后」（`count\tclass`）。
- 所有摘要行统一加 `#` 前缀标记（整段作为注释块，便于日志解析时跳过）。
- 重新编译 release 版，字体修复与本次格式调整均已在 release 生效。

## 2026-09-24（三续）

### 修复染色体图标签 + 新增 NLR 类型图

- 修复 02-chromosome-nlrs.png 横轴染色体名显示异常：`SegmentedCoord::key_points()` 生成的是
  `SegmentValue::CenterOf`，原 formatter 只匹配 `Exact` 导致标签全为空；现改为同时处理
  `Exact`/`CenterOf`/`Last`，并对标签数量加 `n.min(60)` 上限避免过多 contig 时重叠。
- 新增 03-nlr-types.png：NLR 结构域类型计数柱状图（按数量降序）。
  - `nlr-report` 的 `RunStats` 增加 `nlr_type_counts` 字段，`collect()` 按 `domain_string` 汇总。
  - `nlr-plot` 新增 `plot_nlr_type_counts()`。
  - `nlr-cli` 的 `write_outputs` 接入第三张图。
- 验证：`cargo test --workspace --locked` 全部通过；release 运行生成三张非空 PNG
  （非白像素比例约 0.19 / 0.13 / 0.17），02、03 横轴标签区域均有文字像素。

## 2026-09-30

### 流式切片 overlap 由 5000 降为 2000 + 全库错误排查

- 改动：`crates/nlr-cli/src/lib.rs` 的 `RunConfig::new`，`overlap: 5000 -> 2000`
  （`fragment_length` 保持 20000，step 15000 -> 18000，减少约 16.7% 的重复扫描碱基）。
- 依据：组装在 **contig 层面**完成（所有 fragment 的 motif 先按 `dna_sequence_id` 汇总再
  `assemble`），overlap 只负责保证每个 motif 至少被某个 fragment 完整翻译。内置库最大 motif
  宽 50 aa = 150 bp，加上六框翻译末端约 5 bp 截断，overlap 下界约 165 bp；2000 bp 留约 13 倍余量。
- 验证（`cargo test --workspace --locked` 44 项全部通过）：
  - 用改前二进制（overlap=5000）与改后二进制在同一基因组上做 A/B。
  - 黄瓜（Cucumis_sativus）：NLR 位点 74 / 74 一致，complete 56 / 56 一致；
    motif 总数 2061 -> 2008；仅 `5_nlr10` 多出一个 motif_1（P-loop），域类别不变。
  - 大豆（Gmax_508_v4.0）：NLR 位点 598 / 598 一致，complete 366 -> 365；
    motif 总数 13003 -> 12694（-2.4%）；16 个位点的 motif 组成变化。
  - 结论：overlap 变更**不是结果中性**的，位点级有约 0.2% 的边界扰动。
- 根因（重要，属既有设计缺陷）：`nlr-cli::scan_one_fragment` 把 signature 预过滤
  （`has_signature`）作用在**单个 fragment 的 motif 串**上，匹配不到就丢弃该 fragment 的全部命中。
  fragment 边界一变，同一段序列落入的 fragment 及其邻接 motif 就变，预过滤结果随之改变
  （实测某 fragment 在新边界下整段被丢弃、motif 归零）。由此结论：overlap 越大 recall 越高，
  原 5000 在一定程度上是在补偿该设计缺陷。
  - 典型例：`Gm16_nlr56` 在旧档含 3 个尾部 motif_11（LRR）为 TIR-NBARC-LRR（complete），
    新档丢失该片后降为 TIR-NBARC（incomplete）。
- 同时排查出的其他问题（详见当次会话报告，均未改动代码）：
  1. `load_checkpoint` 用 `AssembleParams::default()` 重建组装参数，未传入 `config.assemble`，
     导致 `--checkpoint` + `--relaxed-seed` 二次运行会丢失 relaxed 播种。
  2. `write_outputs` 用 `nlr_seq::fasta::read_all` 一次性读入整条基因组做位点抽提，峰值内存
     与基因组同阶（大豆 950 Mb 基因组实测峰值 948 MB），与 README 的「内存与基因组大小无关」矛盾。
  3. `.motifs.bed` 行序不可复现：`all_motifs` 遍历 HashMap，同一二进制同一输入两次运行
     内容相同但行序不同（`.nlr.txt` 确定）。
  4. GFF 头写 `##gff-version 2`，而 README 表述为 GFF3，二者不一致。
  5. `.stats.tsv` 表头为 3 列（`#section/key/value`），chromosome 行为 5 列；
     `nlr_type_counts` 已统计但从未写入 TSV。
  6. `cargo clippy --workspace --all-targets -- -D warnings`（AGENTS 中列出的命令）会因约 13 处
     既有 lint 失败；CI（rust.yml）只跑 build+test，故当前 CI 仍为绿。
  7. `##date` 用 UTC 而非本地系统时间（Asia/Shanghai 会差 8 小时）。
  8. 文档漂移：`nlr-cli/src/main.rs` 顶部注释仍描述已删除的 `-g/-b/-p/--stats` 等参数；
     `nlr-plot` 顶部注释仍写「两张图」；`RGB_COLORS` 中 21-28 的注释把 TIR 类（25-28）
     写成了 RNL CC。

### 修复 signature 预过滤的切片边界依赖（保留 overlap=2000）

- 目标：保留 overlap=2000 的提速，同时让结果不再依赖切片边界，且不比改前掉 recall。
- 改动（`crates/nlr-cli/src/lib.rs`）：
  - `scan_one_fragment` 不再做预过滤，只负责扫描 + 坐标映射（去掉 `signature_def` 参数）。
  - 新增 `filter_motifs_by_signature`：在排序 + 去重之后，按 (contig, strand, frame) 分组，
    再按相邻 motif 间距 `cluster_gap` 切成局部簇；簇内命中任一 signature（或满足柔性播种）
    才保留整簇。簇边界由 motif 自身坐标决定，与切片相位无关。
  - `cluster_gap` 取 `RunConfig::fragment_length`（20 kb），与原单个 fragment 的判定尺度对齐。
- 关键取舍（实测三档）：
  - contig 级（几乎不过滤）：黄瓜 motif 2061 -> 25616（12 倍），位点 74 -> 78，
    噪声命中大量进入组装，判定为不可用而放弃。
  - 局部簇、`cluster_gap` = seed 合并距离 10 kb：黄瓜 74 位点 / complete 56 与改前完全一致；
    大豆 complete 366 -> 365（`Gm16_nlr56` 的 3 个 motif_11 与主基因体**不同读码框**，
    距最近同框 signature 簇 15.9 kb，被断开丢弃），存在边际 recall 损失。
  - 局部簇、`cluster_gap` = fragment_length 20 kb（最终采用）：黄瓜 74 / complete 57，
    大豆 600 位点 / complete 367（改前 598 / 366），`Gm16_nlr56` 恢复为 TIR-NBARC-LRR。
- 验证：
  - 边界无关性：同一源码分别以 overlap=5000 与 overlap=2000 编译，在黄瓜上
    `.nlr.txt` / `.nlr.bed` / `.stats.tsv` / `.loci.fasta` / `.nbarc.fasta` **逐字节一致**，
    GFF 仅 `##date` 不同，`motifs.bed` 排序后内容一致。
  - 回归：`cargo test --workspace --locked` 全部通过（29 个测试二进制，无 failure）；
    clippy 未新增告警（既有 13 处告警未触碰）。
  - 过滤强度仍保留：黄瓜 motif 25616（contig 级）-> 2403（约 9.4% 留存），与改前 2061（8.0%）同量级。
- 结果：overlap 2000 的约 16.7% 扫描量下降得以保留，同时输出与切片边界解耦、可复现。

## 2026-10-01

### 发布 v1.3.0（推送 + GitHub Release）

- 现状核对：`Cargo.toml` / `Cargo.lock` / `README.md`（badge 与 GFF 示例）/ `docs/index.html` /
  `conda/fastnlr/meta.yaml` 均为 1.3.0，但远端 tag 只有 v1.0.0~v1.2.0，即 1.3.0 尚未发布；
  因此本次直接发布 v1.3.0，不需要改版本号。
- 提交：
  - `365e0d8 fix: make signature pre-filter independent of fragment boundaries`
    （含 overlap 5000 -> 2000 与预过滤粒度改造，见上一条）
  - `98aed94 build: update conda source hash for v1.3.0`
- 标签与发布：创建注解标签 `v1.3.0`（信息 `FastNLR v1.3.0`，指向 365e0d8）并推送，
  `.github/workflows/release.yml` 自动建 Release 并上传三平台二进制与 sha256：
  `fastnlr-x86_64-unknown-linux-gnu.tar.gz` / `fastnlr-x86_64-pc-windows-msvc.zip` /
  `fastnlr-x86_64-apple-darwin.tar.gz`（Release 页：v1.3.0，已置为 Latest）。
- conda 配方哈希：原 `sha256` 实为 v1.2.0 的 tarball 哈希（版本号却已是 1.3.0）。
  已用 `git archive --format=tar.gz --prefix=FastNLR-1.3.0/ v1.3.0` 计算，并与
  `codeload.github.com` 实际 tarball 交叉验证一致：
  `9ffc1e97778e19aa758704ab4abb39db326c7d47c449821f719dd7827254521d`。
  该哈希需在打标签之后才能确定（tarball 内含 meta.yaml，存在不动点问题），
  故沿用仓库既有做法：标签后再补一个提交修正哈希。
- 验证：CI（rust.yml）两次 push 均 success；Release 工作流 success；
  `gh release view v1.3.0` 确认 6 个资产（3 二进制 + 3 校验和）齐备。

## 2026-10-06

### 更新 bioconda PR #69426 至 v1.3.0

- 背景：PR bioconda-recipes#69426 自 2026-09-21 创建后自动检查全绿
  （Lint / Linux Tests / OSX-64 Tests / build and test (ARM) / Summary 均 success），
  但两周多无人评审，且配方仍停留在 1.1.1。
- 改动（commit 76e2169，分支 CropCoder/bioconda-recipes:add-fastnlr-1.1.1）：
  - version 1.1.1 -> 1.3.0，sha256 更新为
    `9ffc1e97778e19aa758704ab4abb39db326c7d47c449821f719dd7827254521d`
    （本次重新独立核验 v1.3.0 tarball，与 10-01 记录一致）。
  - 删除空的 run: 段与 host: pkg-config。
  - build.sh 的 cargo install 增加 --bins。
  - test.commands 增加 `fastnlr --help`。
- 修正（commit 5536b20）：首轮推送时一并删掉了 run_exports，触发 bioconda lint
  的 missing_run_exports（ERROR），导致 Linux/OSX 构建被 skip。已恢复
  `run_exports: - {{ pin_subpackage("fastnlr", max_pin="x") }}`，重推后 lint 通过。
  结论：bioconda 的 missing_run_exports 是对所有 recipe 的强制要求，CLI 工具
  （无下游包）也不能省略，max_pin 按语义化版本取 "x"。
- 同步：仓库内 conda/fastnlr/{meta.yaml,build.sh} 做相同清理，保持与 PR 一致。
- 待办：PR 标题/正文仍显示 1.1.1，需在网页端改为 1.3.0（gh 所用 token 无
  bioconda/bioconda-recipes 写权限，updatePullRequest 返回 403）；催评审评论已由作者
  于 2026-10-06 发出，后续等 reviewer 接手，CI 全绿前不要再推新 commit。

## 2026-10-07

### 补齐 crate 文档 + crates.io 发布 1.3.0 + 收敛为仅 bioconda

- 新增 9 份 crate 级 README（`crates/*/README.md`，共约 1050 行），统一结构：定位 /
  工作区位置 / 模块与算法 / 公开 API / 示例 / 保证 / 测试 / 许可。内容按源码实际签名逐项
  核对，并修正了模块注释与实际不符之处：seed 组合实为 13 条（注释写 11）、`is_seed` 对
  含外部库 motif 的组合改用连续子串匹配、`RGB_COLORS` 中 25-28 的注释误标为 RNL CC
  （实为 TIR）。
- Cargo 打包元数据：7 个 crate 的内部 path 依赖版本由 `1.0.0` 改为与工作区一致的
  `1.3.0`；旧值会让发布到 crates.io 的包声明过宽的依赖范围。README 由 Cargo 自动识别，
  无需改 `Cargo.toml` 的 `readme` 字段（`cargo package --list` 已确认包含）。
- conda：`conda/README.md` 移除 conda-forge 提交流程，明确 bioconda 为唯一目标频道；
  `conda/fastnlr/{meta.yaml,build.sh}` 与 bioconda PR（5536b20）逐字节一致。
- 提交并推送（origin/main `98aed94..36293b0`）：
  - `ab2b592 docs: add per-crate README documentation`
  - `5e399df build: sync conda recipe with the bioconda submission`
  - `36293b0 build: align internal dependency requirements with the workspace version`
- crates.io 发布 1.3.0：此前仅 7 个库 crate 停留在 1.1.1，`nlr-plot` 与 `fastnlr` 从未发布；
  本次按依赖自下而上依次发布 9 个包（nlr-core / nlr-config / nlr-seq / nlr-scan /
  nlr-assemble / nlr-output / nlr-report / nlr-plot / fastnlr）。
  验证：`cargo publish --dry-run -p nlr-core` 通过；发布后 `cargo install fastnlr --version
  1.3.0` 成功安装并输出 `fastnlr 1.3.0`（说明内嵌 motif 数据已随包分发）。
- bioconda PR #69426：作者已合并 master（head `fd9080ae`），Lint / Linux Tests /
  OSX-64 Tests / build and test (ARM) / Summary 全部 success，等待 reviewer 接手；
  PR 标题仍显示 1.1.1，需在网页端改为 1.3.0。
- 按用户要求不再提交 conda-forge，仓库内仅保留"不提交 conda-forge"的说明。

## 2026-10-08

### 修复 BUG-1（位点内 motif 重复）与 BUG-2（NB-ARC 比对错位）

- 背景：核查 1.3.0 时发现 `merge_seeds` 合并重叠种子时会把同一 motif 重复并入位点
  （`MotifList::add_motifs` 只 extend+sort 不去重，`remove_redundant_motifs` 是死代码），
  进而导致 BED12 块重复/重叠、`.nlr.txt` motif 列表虚高；同时
  `write_nbarc_alignment_fasta` 的游标只检查当前位置，遇到重复 motif 后**后续所有
  NB-ARC 块全部退化为 gap**。
- 改动：
  - `crates/nlr-assemble/src/lib.rs::assemble()` 命名循环内新增
    `list.remove_redundant_motifs()`（组装结果唯一出口，按 `(id, dna_start, dna_end)`
    保留排序后首条）。
  - `crates/nlr-output/src/lib.rs::write_nbarc_alignment_fasta()` 比对游标由"只看当前位置"
    改为"从游标起向后查找期望 motif"，命中则推进、缺失才补 gap；干净位点输出等价。
  - 新增 5 个回归测试：`nlr-assemble` 重叠种子去重、`nlr-core` 去重保留首条、
    `nlr-output` 重复 motif 下比对不退化 + 缺块仍正确补 gap、`nlr-cli/tests/outputs.rs`
    真实二进制端到端（共识序列合成 495 bp 基因组，断言 BED 块升序不重叠、txt 列表长度与块数一致）。
- 验证：
  - `cargo test --workspace --locked` 47 项全部通过（原 42 + 新 5）。
  - 反向验证：同样的比对输入在 crates.io 未修复的 1.3.0 上为"长度 165 / 非 gap 115"，
    新测试断言 gap=0 必然失败，确认是有效回归。
  - 真实数据（玉米 chr1 300 kb 切片，同一输入 A/B）：
    motif 条目 17 -> 14、BED 重叠块 3 组 -> 0、NB-ARC 比对非 gap 50 -> 165。
  - `cargo clippy --workspace --all-targets` 仍为 23 条既有告警，未新增。
- 提交（本地，未推送）：`93f436f fix: deduplicate motifs in assembled NLR loci`、
  `4ba22d1 fix: make NB-ARC alignment robust to extra motifs`。
- 未处理：BUG-3（mot/store 最大 id 不一致 panic）、BUG-4（CDF 覆盖不足产生假阳性）、
  BUG-5（dedup_adjacent 三段重复残留）、BUG-6（非 Linux 峰值内存为 0）以及 checkpoint
  参数丢失、read_all 内存、GFF 版本、stats 表头、UTC 时间等既有问题；待全部修完后统一发 1.3.1。

### 发布 v1.3.1（含 BUG-1 / BUG-2 修复）

- 版本号统一更新：`Cargo.toml`（workspace 1.3.0 -> 1.3.1）、7 个 crate 的内部 path 依赖
  版本、CLI banner 与 `--version` 长版本、GFF `##source-version FastNLR V1.3.1`、
  `README.md`（badge 与 GFF 示例）、`docs/index.html`（页头版本与下载链接）、
  `conda/fastnlr/meta.yaml`；`Cargo.lock` 用 `cargo update --workspace` 同步 9 个包版本。
- 提交与推送：
  - `f4bbc94 build: bump version to 1.3.1`
  - `56ae455 build: update conda source hash for v1.3.1`
- 标签与 Release：创建注解标签 `v1.3.1`（指向 f4bbc94）并推送，`.github/workflows/release.yml`
  自动建 Release 并上传 6 个资产（3 平台二进制 + 3 个 sha256）。
- conda 配方哈希：v1.3.1 tarball sha256 =
  `1d47a7148950f821b94a9da357fb6d58086968e13c8f34bdf03aef479bba7a3a`
  （`/archive/refs/tags/` 与 `codeload.github.com` 两路下载交叉核验一致，大小 1,387,740 B）。
  同 v1.3.0，tarball 内含 meta.yaml 存在不动点问题，故标签后再补提交修正哈希。
- 验证：`cargo test --workspace --locked` 47 项全通过；CI（rust.yml）main 分支 success；
  Release 工作流 success；下载 Linux 资产 `sha256sum -c` 通过且 `fastnlr --version` 输出 1.3.1。
- 未做（按本次请求范围）：crates.io 1.3.1 尚未发布（9 个包仍为 1.3.0），bioconda PR #69426
  仍为 1.3.0，待确认后再同步。

### 发布 v1.3.2（修复 BUG-3~6）

- BUG-3（profile 不一致 panic）：`nlr-config::load_from_str` 增加一致性校验，mot.txt 声明的
  motif 若在 store.txt 没有 CDF 表则返回 `InvalidData` 错误并列出缺失 id；CLI 实测
  `error: run failed: store.txt has no CDF table for motif_2 ... `，exit code 2，不再 panic。
- BUG-4（CDF 覆盖不足假阳性）：两处收紧——
  - `cdf()` 越界查询改为返回表内最后一个（最小）p 值（真实 p 的保守上界），不再返回 0.0；
  - `cdf_from()` 对 store.txt 中缺行的分数段默认填 1.0（不显著）而不是 0.0；
  - `score_thresholds()` 在表内无任何分数达标时返回 `i32::MAX`（没有窗口能通过），
    替换原来的 `arr.len()`。
  内置库经核对 0..max 无缺行、表尾 0.0 均为显式写入，故行为中性。
- BUG-5（dedup_adjacent 三段重复残留）：改为删除后不推进索引，连续 3 条以上相同命中全部收敛为 1 条。
- BUG-6（非 Linux 峰值内存为 0）：`read_peak_rss_mb()` 返回 `Option<f64>`，Linux 走
  `/proc/self/status`，macOS 走 `libc::getrusage`（macOS 下新增目标平台依赖），其它平台返回
  `None`，摘要显示 `n/a` 而非误导性的 `0.00 MB`。
- 新增 6 个测试（config 3、dedup 1、摘要格式与 Linux 峰值 2），`cargo test --workspace --locked`
  53 项全通过。
- 回归验证：
  - 内置库结果中性：同一 300 kb 玉米切片在修复前后 `.nlr.txt/.nlr.bed/.motifs.bed/
    .loci.fasta/.nbarc.fasta/.stats.tsv` 逐字节一致。
  - 平台分支：`cargo check --target x86_64-apple-darwin` 与
    `--target x86_64-pc-windows-msvc` 均通过；Release 工作流实际构建出 macOS 资产。
  - clippy 23 -> 22 条（顺手清掉 1 条既有告警，未引入新告警）。
- 提交：`09fe71c`（BUG-3/4）、`26fc012`（BUG-5）、`5a13241`（BUG-6）、
  `5048551`（bump 1.3.2）、`65d5859`（conda 哈希）。
- 发布：注解标签 `v1.3.2`，Release 6 个资产齐备；Linux 资产 sha256 校验通过且输出 1.3.2；
  源码 tarball sha256 = `7909bcebd36fccde5a949c22b473fd528f5fc2a7d182f4bfdd90a8c16ee824a4`
  （archive 与 codeload 双路一致，1,392,242 B）。
- 仍待处理（未包含在 1.3.2）：checkpoint 组装参数不回放、位点抽提 `read_all` 内存与基因组同阶、
  `.motifs.bed` 行序不可复现、GFF `##gff-version 2`、stats 表头列宽与 `nlr_type_counts` 未落盘、
  `##date` 用 UTC、文档漂移；crates.io 与 bioconda 仍停留在 1.3.0。

### 分发渠道同步到 1.3.2（crates.io + bioconda）

- crates.io：按依赖顺序发布 9 个包到 1.3.2（nlr-core -> nlr-config / nlr-seq -> nlr-scan /
  nlr-assemble / nlr-output / nlr-report -> nlr-plot -> fastnlr）。发布前对 nlr-core 做
  `cargo publish --dry-run` 校验；发布后用
  `cargo install fastnlr --version 1.3.2` 实测安装成功并输出 `fastnlr 1.3.2`。
  crates.io API 核对 9 个包 latest 均为 1.3.2（1.3.1 未发布，直接由 1.3.0 跳到 1.3.2）。
- bioconda PR #69426：分支 `CropCoder/bioconda-recipes:add-fastnlr-1.1.1` 从 `fd9080a`
  更新到 `b4aeb4c Update fastnlr to 1.3.2`，配方改为 version 1.3.2 + sha256
  `7909bceb...ee824a4`，与仓库内 `conda/fastnlr/` 逐字节一致；build number 重置为 0。
  验证：Lint success，Linux Tests / build and test (ARM) 已启动，等待整体跑绿后由 reviewer 接手。
- 遗留：PR 标题仍为 "Add fastnlr 1.1.1"，需在网页端改为 1.3.2（当前 gh token 对
  bioconda/bioconda-recipes 无写权限）。

## 2026-10-09

### 新功能：蛋白序列模式（HMMER 结构域扫描）并发布 v1.4.0

- 需求：支持蛋白 FASTA 输入，用 hmmer-pure-rs 做结构域扫描识别 NLR 特征结构域并给出结果。
- 方案选型（先做可行性验证，全部在 /tmp 完成）：
  - `hmmer-pure-rs 0.7.4`：纯 Rust 的 HMMER 3.4 移植，仓库内 LICENSE 为 BSD-3-Clause（与
    GPL-3.0 兼容，crates.io 元数据未声明 license 字段）；本机编译通过，release 下约
    1300 次 模型×序列 扫描/秒。
  - 精选 4 个 Pfam 模型并核对 accession：**PF00931**(NB-ARC, GA 23.5/23.5)、
    **PF01582**(TIR, 21.3/21.3)、**PF18052**(Rx_N/CNL 型 CC)、**PF05659**(RPW8/RNL 型 CC)。
    注意 PF01560 在当前 Pfam 已不是 TIR（是 HCV_NS1），不可沿用旧经验编号。
  - **LRR 不走 HMM**：7 个 Pfam LRR 模型联合在 RefPlantNLR 上灵敏度仅 47.5%（假阳性 0.2%），
    而 FastNLR 现有 PWM motif（motif_9/11）在同集上为 99.5%，故采用 HMM(核心域) + PWM(LRR) 混合。
- 实现：
  - 新 crate `nlr-domain`（第 3 层）：内置 4 个 HMM（`crates/nlr-domain/data/`，共 277 KB，
    include_str! 内嵌，附 README 与 `fetch_hmms.py` 复现脚本）；`HmmLibrary`（同 accession
    覆盖、新 accession 追加）、线程内 `Scanner`（OProfile 非 Clone/Send，故用 rayon map_init
    按线程构造）、`classify()` 架构判定（nterm/class/complete/is_nlr）。
  - CLI：新增 `--protein` 与 `--hmm`；蛋白模式与 `--flank/--relaxed-seed/--flexible-seed/
    --require-ploop/--signature/--seed-combination/--checkpoint` 互斥并给出明确报错。
  - 输出：`{prefix}.domains.tsv`（HMM + PWM 逐条命中，含 GA 通过标记）、`{prefix}.nlr.tsv`
    （候选位点架构）、`{prefix}.summary.txt`、`{prefix}.plots/{01-nlr-types,02-domain-counts}.png`；
    不生成 GFF/BED/loci/nbarc。
  - `nlr-plot` 抽出通用 `plot_counts()` 供两种模式复用。
- 测试：新增 7 个（nlr-domain 4 + CLI 端到端 3），`cargo test --workspace --locked` 60 项全通过；
  含 golden bitscore（±0.5 bit）与打乱序列阴性对照。
- 验收回归（RefPlantNLR 400 条，真实二进制）：NB-ARC 398/400(99.5%)、TIR 89/91(97.8%)、
  CC 246/263(93.5%)、LRR 391/394(99.2%)，打乱对照假阳性 0/400；400 条蛋白 8 线程 0.465 s。
- 打包风险核实：`cargo-bundle-licenses`（现装 v4.2.0）实测 exit 0，正确识别 hmmer-pure-rs 为
  BSD-3-Clause 并抓取其 LICENSE 文本，conda/bioconda 构建不受影响。
- 发布 v1.4.0：
  - 提交 `3e3229e feat: add protein input mode with HMMER domain scanning`、
    `dbc6ed8 build: bump version to 1.4.0`、`467cb92 build: update conda source hash for v1.4.0`。
  - 注解标签 `v1.4.0` → Release 6 个资产齐备；tarball sha256 =
    `e12d89da56912c92e11aaedef527da5304fa066e5a4c75f9a1bcbf0196f3c16d`（双路一致）。
  - crates.io：**10 个包**全部发布 1.4.0（含首次发布的 `nlr-domain`）；`cargo install
    fastnlr --version 1.4.0` 实测可装，且安装版 `--protein` 在三条约 1 kb 的 NLR 夹具上输出
    2 条 CC-NBARC-LRR + 1 条 TIR-NBARC-LRR，说明 HMM 数据已随包分发。
  - bioconda PR #69426：分支推进到 `8e209a8 Update fastnlr to 1.4.0`，配方与仓库内副本逐字节一致；
    Lint success，Linux/ARM 构建进行中。PR 标题仍需在网页端更新（token 无该仓库写权限）。

## 2026-10-09（续）

### 重建 benchmark 测评体系（v1.4.0）并完成 Part 1 准确度评测

- 背景：旧 6 个脚本使用已废弃的 `-o out.txt -p` 接口（实测 `unexpected argument '-p'`），
  09-21/22 的旧数据全部作废；本次按新方案重建协议并重跑。
- 协议文档：新增 `benchmark/README.md`（数据来源、真值映射、测量口径、机时预算、复跑方式）；
  `Script3/4` 标注为 SUPERSEDED（功能被 Script1 的 D2 覆盖）。
- 脚本改造（原地升级为新 CLI）：
  - `Script5`：RefPlantNLR 400 条 CDS 跑四组配置（dna_default / dna_flexible `--flexible-seed 3` /
    dna_extended 28 motif / protein `--protein`），CDS 由 R 按标准遗传密码翻译成蛋白供蛋白模式使用；
    真值 `Structure` 按位映射并折叠相邻重复；输出逐基因表、按类别、按生物学类型、按 20 类原样、
    总体指标共 5 个 CSV。
  - `Script6`：FigA1 分类型召回、FigA2 架构一致率热图、FigA3 N 端混淆矩阵、FigA4 LRR 与完整性。
  - `Script1`：D1 普查（18 物种 × t8 × 3 次）+ D2 线程扩展（5 物种 × {4..256} × 2 次）
    + D3 参数档（3 物种 × 3 档）；`/usr/bin/time -v` 采集 wall/RSS/CPU/IO，并与 `.summary.txt`
    自报值交叉校验；运行前预热页缓存；增量写 CSV 支持断点续跑（`--collect-only`、`--rerun`、
    `--designs=`、`--only-species=`）。
  - `Script2`：FigB1-B8（时间/内存 vs 基因组大小、线程扩展、加速比与效率、吞吐量、
    各物种 NLR 数、参数效应、NLR 类型构成）。
- **Part 1 结果**（RefPlantNLR 400 条，四组配置）：

  | 配置 | 召回 | 架构一致率 | LRR 灵敏度 | 过度切分 |
  |---|---|---|---|---|
  | dna_default | 393/400 (98.25%) | 91.6% | 99.7% | 3 |
  | dna_flexible | 398/400 (99.50%) | 89.4% | 99.7% | 3 |
  | dna_extended | 398/400 (99.50%) | 82.7% | 99.7% | 28 |
  | protein (HMM) | 398/400 (99.50%) | **94.2%** | 99.7% | n/a |

  按生物学类型召回：TNL 86 条四组均 100%；CNL 263 条 99.2-99.6%；NL 34 条 97.1-100%；
  **RNL 8 条从 default 3/8 提升到 flexible / extended / protein 的 8/8**；TN 5 条全 100%。
  未检出：default 7 条（ADR1、ADR1-L1/L2/L3、AtNRG1.1、NRG1、Wed），
  flexible/extended 2 条（Pb1、Wed / ADR1-L3、Wed），protein 2 条（CW9、Pb1）。
  亮点：蛋白模式用 HMM 直接判结构域，架构一致率最高（94.2%），且 RNL 全部召回。
- **Part 2 已启动**（后台 setsid 运行，跨会话存活）：133 个场景（survey 54 / threads 70 / params 9），
  已完成 3 个（冒烟测试的黄瓜 t8 × 3，53.5 s / 206-210 MB / NLR 74 / complete 57 / motif 2403，
  三次完全一致）。1.4.0 比旧版快约 2.4 倍（黄瓜 t8：130 s → 53.5 s），预计总机时约 7 h。
  运行日志 `benchmark/logs/campaign.log`，逐次指标增量写入 `benchmark/metrics/genome_benchmark_raw.csv`。

## 2026-10-10

### Part 2 性能测评收尾（v1.4.0，18 基因组 × 线程矩阵）

- 执行状态：133 个场景全部执行，132 个完成、0 个失败；最后一个
  `Zanthoxylum_bungeanum_t8_sensitive_r1`（4.7 Gb + sensitive）按用户决定在运行约 3 小时后终止，
  raw CSV 记为 `exit_status = 15`（SIGTERM）、指标为 NA，manifest 标注为主动跳过，
  README 中注明 sensitive 档仅在黄瓜（197 Mb）与大豆（991 Mb）上做了表征。
- 记账问题与修正：run_id `{物种}_t{线程}_{模式}_r{重复}` 在 t8+default 档被 D1/D2/D3 三个设计共用，
  CSV 按 run_id 去重导致 13 行被覆盖（5 个线程物种的普查 t8 行、3 个参数物种的 t8 default 行）；
  测量未丢失。`Script2` 已改为按“配置”（`threads==8 && param_mode==default`）归并普查点与线程曲线，
  每个物种恢复 3–4 次重复；协议文档 `benchmark/README.md` 增补“已知记账问题”与“执行状态”两节。
- **性能结论**（全部为串行 warm-cache 实测，`/usr/bin/time -v`，摘要自报值交叉校验通过）：
  - 普查（18 物种 × t8 × 3 次，中位数）：197 Mb 黄瓜 53.5 s / 210 MB → 4.7 Gb 花椒 1229.9 s / 4625 MB；
    单次普查总机时约 1.67 h。**峰值内存与基因组大小基本 1:1（比值 0.97–1.15）**，证实位点抽提
    `read_all` 的内存问题，README 中“内存与基因组规模无关”的表述需要改写。
  - 线程扩展（自 t4 起）：黄瓜 15.4×（121→8 s）、马铃薯 13.3×、大豆 13.9×、花椒 18.9×（2451→130 s）；
    **花生在 256 线程崩塌**（t128 89 s → t256 459 s，两次重复均复现，CPU 占用反而更高），
    属高线程数资源争用。
  - 参数档：relaxed 档位点略降（黄瓜 74→67、大豆 600→513、花椒 616→575）；
    **sensitive 档病态放大**（黄瓜 53.5 s→2189 s、位点 74→8833；大豆 270 s→9079 s、位点 600→48241，
    RSS 968 MB→5178 MB），不可作为默认。
- 产物：`benchmark/figures/FigB1–B8.pdf`（Nature 风格、Arial、cairo_pdf）；
  `benchmark/metrics/genome_benchmark_summary.csv`（每物种中位数）、`genome_benchmark_raw.csv`（120 行）、
  `run_manifest.csv`；修正了两处绘图 bug（B8 调色板不足、B3–B5 物种配色未匹配标签）。
- 准确性部分（Part 1）结论见上一节：四组配置召回 98.25%/99.50%/99.50%/99.50%，
  蛋白模式架构一致率最高（94.2%），RNL 由默认 3/8 提升到三种现代方案的 8/8。
- 遗留（下次整理）：`benchmark/accuracy/protein_dna_mode_misrun/` 为脚本 bug 造成的误跑残留（约 40 KB），
  待用户确认后清理；benchmark 目录仍在 `.gitignore` 内（本地资产）。

### FigA1 改版与 "Other" 类别溯源（准确性类型分组修正）

- 改版需求：纵轴由召回百分比改为**检出基因数**、去掉百分比刻度、改用项目主题色的浅色低对比配色、
  缩小尺寸并调整比例（180×100 mm → 120×72 mm）、图例改为底部简短说明
  （Default (20 motifs) / Flexible seed / Extended (28 motifs) / Protein (HMM)），
  并在标题下加一行提示 “Bar height = genes detected; n = genes in that class”。
- 溯源发现：原图 "Other (n=4)" 并非生物学类别，而是我分组规则的兜底桶——这 4 条
  （NbPRFa、NbPRFb、Prf、R1）的 Structure 以 **B** 开头（BCCNL / BCNL），而规则只识别 C/T/R/N/O。
  查 xlsx 的 Key 页得到权威定义：**B = PF12061（Late blight resistance protein R1 型 N 端域）**，
  C = PF18052(Rx_N/CC)、T = PF01582(TIR)、R = PF05659(RPW8)、N = SSF52540(NB-ARC)、
  L = LRR、O = 其他整合域；这 4 条的 Subclass 注明为 `NLR;CNL;NRC`，属 NRC 支系 CNL。
- 修正：`Script5` 的生物学类型规则把首字母 B 与 C 同等对待（`^[CB]` → CNL），"Other" 桶消失，
  CNL 由 263 增至 **267**；对应的检出数（Default/Flexible/Extended/Protein）为
  266/265/266/265，其余类型不变（TNL 86/86、RNL 3/8/8/8、NL 33/34/33/34、TN 5/5/5/5）。
  总体指标不变：召回 98.25%/99.50%/99.50%/99.50%，架构一致率 91.6%/89.4%/82.7%/**94.2%**。
- 附带发现（值得进 1.5.0）：蛋白模式对 4 条 B 型 CNL 中的 3 条（NbPRFa、NbPRFb、R1）判为
  `NBARC-LRR`，未识别其 N 端结构域；DNA 模式 4/4 正确。原因是 `nlr-domain` 目前只内置
  PF18052(Rx_N) 与 PF05659(RPW8) 两个 CC 模型，缺 **PF12061**。建议下一版把 PF12061 加入内置
  HMM 集合（类别 CC），可补上 NRC 支系的 N 端结构域判定。

## 2026-10-11

### 内置 HMM 集合加入 PF12061（R1 型 N 端结构域，类别 CC）

- 背景：上一条发现蛋白模式在 NRC 支系 CNL（金标准 `Structure` 以 B 开头，B = PF12061）上漏判
  N 端结构域。本次把 PF12061 加入内置模型集合。
- 模型：`PF12061.14`，Pfam 短名 *NB-LRR*，描述 “Late blight resistance protein R1”，
  长度 295，GA 28.7/28.7，HMM 文本 138 KB；来源 InterPro REST API（同既有四个模型）。
- 代码改动：
  - `crates/nlr-domain/data/PF12061.hmm` + `data/README.md`（表格与选型理由）+ `fetch_hmms.py`
    （MODELS 增加 PF12061）。
  - `EMBEDDED_MODELS` 由 4 个增至 5 个（名称 `R1-type`，类别 `Cc`）；
    `category_for_accession` 增加 `PF12061 => Cc`；CLI 启动日志的模型列表同步。
  - 测试夹具新增 Prf（UniProt Q96485，番茄，1824 aa，与金标准条目一致）；
    `tests/scan.rs` 断言模型数 5、Prf 的 PF12061 命中 161–267 位、bitscore 34.0（±0.5）；
    `crates/nlr-cli/tests/protein_mode.rs` 夹具数 3→4 并新增 Prf 的类别断言。
  - `cargo test --workspace --locked` 60 项全部通过。
- 效果（重跑 `Script5 --rerun` 与 `Script6`）：
  - 召回不变：98.25% / 99.50% / 99.50% / 99.50%（PF12061 只影响架构判定，不影响检出）。
  - **蛋白模式架构一致率 94.2% → 95.0%**；4 条 B 型 CNL（NbPRFa、NbPRFb、Prf、R1）
    现在全部判为 `CC-NBARC-LRR`（此前 3 条被判成 `NBARC-LRR`）。
  - 蛋白模式剩余架构不一致 20 条：CC-NBARC-LRR→NBARC-LRR 15 条、TIR-NBARC-LRR→NBARC-LRR 2 条、
    多域串不一致 2 条、CC-NBARC-LRR→CC-NBARC（丢 LRR）1 条。
- 待办：本次改动尚未发版；若要分发（crates.io / bioconda / 二进制）需升 1.4.1 并同步配方哈希。

### bioconda 上架后的文档同步（README 与官网主页）

- 背景：FastNLR 1.4.0 通过 bioconda 审核并上架
  （`https://anaconda.org/bioconda/fastnlr`），安装方式推荐 `conda install bioconda::fastnlr`。
  核对 Anaconda API：`latest_version = 1.4.0`，四个子平台 linux-64 / linux-aarch64 / osx-64 / osx-arm64 均有产物。
- README 更新：
  - 徽章区新增 bioconda 版本、下载量、平台三个徽章（`anaconda.org/bioconda/fastnlr/badges/*.svg`）。
  - 新增 `## Install` 章节，把 bioconda 列为推荐安装方式（含 conda-forge 通道顺序与
    `bioconda::` 前缀说明、mamba 等价命令、包页与 `conda/` 配方位置），原二进制/源码安装顺延为 B/C 选项。
  - 功能列表新增“一行安装”条目；蛋白模式模型列表补齐 R1-type。
  - 修正陈旧数字：测试数 40 → 60（实测 `cargo test --workspace --locked` 通过 60 项）。
  - 重写性能章节为 18 基因组实测数据（197 Mb 黄瓜 53.5 s/210 MB → 4.7 Gb 花椒 1230 s/4625 MB，
    吞吐 3.7 Mb/s，峰值内存与基因组约 1:1），并新增 RefPlantNLR 召回与蛋白模式架构一致率小结；
    删除与实测矛盾的“内存与基因组大小无关”表述。
- 官网 `docs/index.html` 更新：首屏加 bioconda 安装条与 Linux/macOS 平台标注；导航与区块
  `Quick Start` 改为 `Install`（新增 Conda 页签，默认激活）；特性网格补“One-line install”卡片；
  蛋白模式卡片补 R1-type 与 `--hmm`；修复卡片与“实现说明”里已废弃的 `-a/-f/-g/-b/-m/-c` 旧短选项；
  性能区统计带改为 18 基因组/60 测试/6 阅读框/4 平台，表格换成实测四物种数据，新增召回表与 95% 架构一致率；
  页脚增加 bioconda 链接。
- 文档结构按要求调整为先简介、再安装与用法、最后特性：README 章节顺序改为
  `Install → Usage → Features → …`（原 `Features` 前置移到 `Usage` 之后，去掉与安装章节重复的
  “一行安装”条目，并在安装章节末尾加用法指引）；官网区块与导航同步改为
  `Install → Usage → Features → Architecture → Performance → Correctness`，并重排相邻区块的
  底色 class 以保持浅/深交替。
- 验证：`docs/index.html` 经 HTML 标签配对解析无未闭合/错配；用无头 Chrome 渲染首屏、特性网格、
  安装区、性能区截图逐一核对版式，无重叠或溢出。

### 代码整理与 v1.4.1 收尾发版准备

- 版本单一来源改造（此前版本号散落在 3 处硬编码）：
  - `crates/nlr-cli/src/main.rs` 的 ASCII banner（`FastNLR v1.4.0`）与 `long_version` 改用
    `env!("CARGO_PKG_VERSION")`，banner 拆为 `concat!` 两段并保持原有排版（含行尾空格）逐字节不变；
  - `crates/nlr-output/src/lib.rs` 的 GFF `##source-version` 改用包版本，输出实测为 `FastNLR V1.4.1`。
- 版本提升到 1.4.1：workspace `version`、16 处 crate 内部依赖约束、`Cargo.lock`（10 个包）、
  README 徽章与 GFF 示例、官网版本徽标与 Release 下载链接、`conda/fastnlr/meta.yaml` 的 `version`。
- 代码质量修复（此前 `cargo fmt --check` 有 23 个文件 86 处 diff、clippy 有约 20 条告警）：
  - `cargo fmt --all` 全量整理，现 `cargo fmt --all -- --check` 通过；
  - clippy 逐条修复：`manual_range_contains`（motif.rs、config 测试）、`should_implement_trait`
    （`DomainCategory::from_str` 更名 `parse`，同步 nlr-cli 调用点）、`explicit_auto_deref`
    （`.map(|v| *v)` → `.copied()` ×2）、`collapsible_str_replace`、
    `needless_range_loop`（SIMD lane 改 `enumerate`）、`unnecessary_cast`、`ptr_arg`
    （`&PathBuf` → `&Path` ×2）、`type_complexity`（测试加类型别名）；
    `arbitrate` 的 `too_many_arguments` 就地保留并加 `#[allow]` 与理由注释。
    现 `cargo clippy --workspace --all-targets -- -D warnings` 零告警。
  - 顺带修正 `crates/nlr-cli/src/main.rs` 顶部已过时的模块文档（旧短选项列表）。
- 未改动项：`benchmark/Script3` 仍用旧 `-p` 参数，但脚本头部已标注 `[SUPERSEDED]`“不再维护、
  保留作历史记录”，故按归档件保留；`benchmark/` 目录在 `.gitignore` 内，不进版本库。
- 验证：`cargo test --workspace --locked` 60 项全通过；`cargo build --release --locked` 成功，
  `fastnlr --version` 输出 `fastnlr 1.4.1`；在 20 kb 合成序列上跑通全流程，GFF 头部为
  `##source-version FastNLR V1.4.1`，9 类输出文件齐备。
- 提交与发布（5 个提交，均为规约式前缀、不含工具署名）：
  `a250eae style:` → `a0fbb41 refactor:` → `20ddd40 build: bump version to 1.4.1`
  → `fe4f1bd docs:` → `0c92263 chore:`，推送到 `origin/main`；随后打注解标签 `v1.4.1` 并推送，
  触发 Release 工作流（run 38116894688，2m48s，成功）。
- Release 资产 6 个齐备（linux-gnu / apple-darwin / windows-msvc 三个 archive + 三个 `.sha256`），
  校验发布版 Linux 二进制 `--version` 输出 `fastnlr 1.4.1`。
- conda 同步：`v1.4.1` tarball sha256 = `ea174c2c3fb8a1e776d29baaae0a864c728978a48271683795905361cf63df25`
  （两次独立下载一致），写入 `conda/fastnlr/meta.yaml` 后以 `13e87a5 build: update conda source hash
  for v1.4.1` 推送。CI（push main）与 GitHub Pages 部署均为 success，线上主页已含
  `conda install bioconda::fastnlr` 与 v1.4.1。
- 遗留（需用户决定，未执行）：① 向 `bioconda/bioconda-recipes` 提 PR 把配方由 1.4.0 升到 1.4.1
  （v1.4.0 走过 PR 流程，本次未动外部仓库）；② crates.io 是否同步发布 1.4.1 的 10 个 crate；
  ③ 是否把 `cargo fmt --check` 与 `cargo clippy -D warnings` 加入 CI（涉及 CI 配置，未擅自修改）。

## 2026-10-12

### 拟南芥 TAIR12 真实基因组评估数据准备与 FastNLR 首次分析

- 数据获取：从 NCBI 下载 TAIR12（`GCA_978657495.1`，Complete Genome）到
  `benchmark/arabidopsis_tair12/raw/`，7 个文件共 74 MB（基因组 37 MB、GFF 14 MB、蛋白 7.1 MB、
  CDS 13 MB、feature_table 2.6 MB + assembly_report + md5checksums），
  MD5 与 NCBI 官方校验和全部一致。未下载 gbff/gtf/gpff/rna（重复且体积大）。
- 数据概况：5 条染色体 142.5 Mb，58,717 个基因（蛋白编码 26,867），39,150 个 mRNA，32,589 条蛋白。
  染色体 accession 映射 OZ408683.1→Chr1 … OZ408687.1→Chr5。
- **注释文件无法直接给出 NLR 清单**（重要结论）：
  - GFF 的 `gene` 只有 `locus_tag/gene_biotype/Name/ID/gene`，`gene_biotype` 无 NLR 类；无 GO 注释；
  - `mRNA` 的 `product=` 是自由文本，且命名体系不一致：RPS4(AT5G45250) 写作
    "putative NAD+ nucleotidase, cyclic ADP-ribose generating"，RPS5(AT1G12290) 写作
    "apoptotic protease-activating factor, helical domain-containing protein"，**不含任何 NLR 关键词**；
  - 常规 NLR 关键词只能命中 54 个基因（漏掉 RPS4/RPS5），改用 TAIR12 自身的
    "Apoptotic protease-activating factor / NAD+ nucleotidase / putative ATP hydrolase" 术语则命中 341 个
    （过检，"putative ATP hydrolase" 一项 204 个，含大量非 NLR ATP 酶）；拟南芥实际约 150 个，
    无阈值可切分。
  - 可用于建金标准的只有：`locus_tag` 里的 AGI 基因号（35,223 个，可精确连接外部策展名单）
    与少量基因符号（`gene=TN2`、`[RESISTANT TO P SYRINGAE 2]`）。
- FastNLR 基因组模式运行（`outputs/fastnlr_dna_default/`）：8 线程 **39.7 s / 162 MB**，
  检出 **207 个 NLR 位点**（完整 139）。类型：TIR-NBARC-LRR 90、TIR 36、CC-NBARC-LRR 30、
  NBARC-LRR 21、TIR-NBARC 19、NBARC 11。染色体分布 Chr1 69 / Chr2 11 / Chr3 19 / Chr4 40 / Chr5 68。
  核对 `is_complete_nlr`（P-loop + LRR）判据得 139，与工具自报一致。
- 可视化：新增 `benchmark/Script7-Plot-TAIR12-NLR-V1.R`（Nature 风格、cairo_pdf、180 mm 宽、
  Liberation Sans 代替未安装的 Arial、图中无中文），输出到 `benchmark/figures/`：
  `Fig11-TAIR12-genome-map.pdf`（染色体分布，正负链分列，四类着色）、
  `Fig12-TAIR12-nlr-composition.pdf`（类型完整度构成 + 染色体重分布）、
  `Fig13-TAIR12-nlr-profile.pdf`（motif 使用谱对数轴 + 位点跨度分布）。
  修复了两处：Fig12b 零值标签显示为 "0"、Fig13a 量纲被 motif_11 压扁（改对数轴）与分组顺序不统一。
- 本机环境核查：Resistify 2.4.0 已装于 conda env `resistify`（蛋白级 NLR 分类器，可作对比工具）；
  HMMER 3.3.2、Java 11、R 4.5.1 可用；bioconda 无 nlr-parser/rgaugury/nlrscape/nlrexpress/drnlr；
  本机无 Pfam-A.hmm。
- 遗留（待用户决定）：金标准路线（策展 AGI 名单 / 结构域推导 / 工具共识）、对比工具范围、
  评估层面（基因组层区间求交 vs 蛋白层直接比对）。`benchmark/` 在 `.gitignore` 内，
  本次脚本与图件均为本地资产，未入库。

### 修复：--flexible-seed 的播种语义由「替换」改为「并集」

- 问题定位（承接金标准漏检分析）：`find_seeds` 的播种判定是 if/else 链，
  开启柔性播种后**只**看 `is_seed_flexible`，精确种子组合被整体顶替；而签名预过滤
  `has_signature_flexible` 本就是 `has_signature || is_seed_flexible` 的**并集**。
  二者不一致导致「motif 簇因命中精确签名被保留、播种阶段却不认精确组合、最终组装不出位点」，
  典型病例就是 RefPlantNLR 的 Pb1（命中精确组合 `[12,2,8]`，但 NB-ARC 类 motif 仅 2 个、
  达不到柔性门槛 3）——默认能检出，加 `--flexible-seed 3` 反而漏检。
- 改动：`crates/nlr-assemble/src/lib.rs` 的 `find_seeds` 一行，改为
  `def.is_seed_flexible(&s, &cfg) || def.is_seed(&s)`，并加注释说明并集语义与理由。
  并集只新增种子、不移除任何原有种子，因此**不传该参数时行为完全不变**。
- 回归测试（`crates/nlr-assemble/tests/assemble.rs` +2 项，共 62 项）：
  ① `flexible_seed_is_additive_not_replacing`——用 Pb1 的形状（精确组合 `[12,2,8]` +
  仅 2 个 NB-ARC motif）断言柔性模式下仍能播种；已实测该测试在未修复版本上**必然失败**
  （`开启 --flexible-seed 不应丢失精确组合能检出的位点`）。
  ② `flexible_seed_still_recovers_rank_order_loci`——反向保护柔性路径不被削弱。
- 文档同步：README 的 `--flexible-seed` 说明与 CLI help 文本由「instead of」改为
  「additionally / additive」；`docs/flexible-seeding-upgrade.md` 补修订说明。
- 验证（修复后二进制，8 线程）：
  - 金标准 RefPlantNLR flexible 档：召回 398/400 → **399/400（99.75%）**，漏检仅剩 Wed；
    过度切分不变（仍是 Prv^1、Prv^2、RPP2A 共 3 个基因）；架构一致率 356/398 → 356/399
    （89.22%，因新召回的 Pb1 架构判别为 CC-TIR-NBARC-LRR、与真值 CC-NBARC-LRR 不符，属预期）。
  - 零副作用确认：RefPlantNLR 的 default 与 extended 档输出、TAIR12 的 default 档输出
    均与修复前**逐字节一致**（extended 不用柔性播种，default 走 else 分支）。
  - TAIR12 真实基因组 flexible 档：位点 176 → **216**；被 flexible 完全漏掉的 default 位点
    **38 → 0**；default 检不到的新位点仍为 6（柔性独有能力未削弱）。
  - 分裂代价：自身重叠/嵌套位点对 default=2、flexible 未修复=1、flexible 已修复=2
    —— 修复后与默认档**持平**，未超出默认档既有的分裂水平；其中 1 对
    （OZ408683.1:29638277-29643710 与 29639806-29640166）与 default 完全相同，属继承而非新增。
  - 性能：flexible 档 39.26 s → 39.36 s（+0.3%，在噪声范围内），峰值内存 162 MB 不变。
  - 质量门禁：`cargo fmt --check` 通过、`cargo clippy --all-targets -- -D warnings` 零告警、
    `cargo test --workspace --locked` 62 项全通过。
- 遗留：`benchmark/metrics/accuracy_summary.csv` 中 flexible 行的召回仍为修复前的 99.50%，
  需要重跑 Script5 才会刷新（未擅自改动既有测评产物）。

### 日志时间戳改为本地时区 + 柔性播种默认开启

- 需求：① 日志前的时间戳改用运行主机时区，格式 `YYYY-MM-DD HH:MM:SS`（秒不保留小数）；
  ② `--flexible-seed` 默认开启，提高召回，同时避免位点分裂与过度切分。
- 时间戳改造：
  - `crates/nlr-cli/src/main.rs` 新增 `LocalTimer`（实现 `tracing_subscriber::fmt::time::FormatTime`），
    `init_logging` 用 `.with_timer(LocalTimer)` 挂载；格式常量 `DATE_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"`。
    此前用的是 tracing 默认 `SystemTime`（UTC、`2026-10-11T06:03:12.917547Z`，带 6 位小数）。
  - `now_date_string()`（GFF `##date` 头）同步改为本地时间，并删除手写的
    `civil_from_days()` UTC 日期换算（约 25 行），统一由 chrono 负责。
  - 新增依赖 `chrono 0.4`（`default-features = false, features = ["clock", "std"]`），
    写入 workspace 依赖表并注明用途；本地 registry 已有缓存，离线可构建。
  - 实测：主机 `Asia/Shanghai`，日志与 GFF `##date` 均为 `2026-10-12 15:18:58`，
    与 UTC（`07:18:58`）相差 8 小时，符合预期。
- 柔性播种默认开启：
  - `--flexible-seed <N>` 保留用于调阈值；新增 `--no-flexible-seed` 回退到仅精确组合。
  - 解析逻辑抽到 `nlr_cli::resolve_flexible_seed()`（含 `DEFAULT_FLEXIBLE_SEED_MIN_NBARC = 3`），
    main 与单测共用同一实现。
  - 启动日志新增一行 `seeding:` 明确当前档位，例如
    `seeding: flexible (min_nbarc=3)` / `seeding: exact seed combinations only`。
  - 校验：`--flexible-seed` 与 `--no-flexible-seed` 互斥；`--require-ploop` 需柔性开启；
    `--relaxed-seed` 与默认柔性互斥（提示改用 `--no-flexible-seed`），避免被静默忽略；
    `--protein` 模式下仍按"基因组专属参数"统一报错（播种校验已用 `!cli.protein` 让位）。
- 新增测试 3 项（`crates/nlr-cli/tests/logging.rs`）+ 单测 1 项（`resolve_flexible_seed`），
  合计 66 项全通过；`cargo fmt --check` 与 `clippy -D warnings` 均干净。
- 验证（金标准 RefPlantNLR 400 条 CDS）：
  - 新默认档召回 **399/400（99.75%）**，旧默认档 393/400（98.25%）；
  - **过度切分不变**（同为 3 个基因：Prv^1、Prv^2、RPP2A）；
  - 架构一致率 91.60% → 89.22%（新召回的基因架构判别不完美，属预期）；
  - `--no-flexible-seed` 的输出与旧默认档**逐字节一致**（可完整回退）；
  - 新默认档与显式 `--flexible-seed 3` **逐字节一致**。
- 验证（TAIR12 真实基因组 142 Mb）：
  - 位点 207 → 216；旧默认档位点被新默认档**完全漏掉 0 个**；
  - 自身重叠/嵌套位点对：旧默认 2、新默认 2 —— **与旧默认持平，未加剧分裂**；
  - 运行 39.3 s / 162 MB，与旧默认档持平；protein 模式输出**逐字节不变**。

### 结构域判别优化 Phase 1：有序域模型（2026-10-12）

- 背景：架构一致率 89.22%（399 检出中 43 例判错）。误差分解——A 类「漏判 N 端域」
  23 例（53%）、B 类「多判 N 端域」9 例（21%，其中 7 例金标准 Subclass 本身写 CNL、
  属参考标注不一致）、C 类「域名串构造伪影」8 例（19%）。Phase 1 只针对 C 类。
- 根因：原 `domain_string()` 仅做「相邻同类合并」，无位置/顺序约束，产出
  `LRR-CC-NBARC-LRR`（单 motif LRR 跑到 CC 之前，4 例）、`CC-NBARC-CC-NBARC-LRR`
  （NB-ARC 之间混入杂散 CC，1 例）等金标准中不存在的形态。
- 改动：
  - `crates/nlr-core/src/signature_def.rs` 新增 `canonical_domain_string(cats, repeat_min_block)`，
    实现「有序域模型」：折叠域块 → 以首个 NB-ARC 为锚点 → N 端域取锚点前最靠近的 CC/TIR →
    锚点起按位置输出，NB-ARC/LRR 全保留，CC/TIR 须同时满足「后随 NB-ARC 块」且
    「块规模 ≥ repeat_min_block」才认定为真实串联重复 → 相邻同类折叠。
  - `DNA_REPEAT_MIN_BLOCK = 2`、`PROTEIN_REPEAT_MIN_BLOCK = 1`：DNA 的一个域块是 PWM 命中
    集合，NB-ARC 之后的**单点** PWM 命中极可能是噪声；蛋白模式的一个块是已过 GA 的 profile
    HMM 命中，本身即域级证据，故不做门槛（否则 RPP2A 的串联 TIR 会被误丢）。
  - `crates/nlr-domain/src/lib.rs` 的 `classify()` 改为调用同一函数，DNA 与蛋白两条通路
    从此共用一套判定口径（此前各写一套，实测同一基因两模式判定不一致率 11.8%）。
- **被数据否决的方案（记录以免重犯）**：原计划给域块设「最小规模门槛」（CC/TIR ≥2、LRR ≥2）。
  实测 CC 域块 249 个中有 85 个只命中 1 个 motif、LRR 块有 34 个只命中 1 个，
  加门槛后金标准一致率从 90.5% 暴跌到 73.7%（CC 门槛）/ 83.7%（LRR 门槛），故明确不加。
- 锚点选择亦有实测依据：RPM1 的域块 `CC(1) NBARC(1) CC(1) NBARC(8) LRR` 中，两个 NB-ARC 块
  被杂散 CC 隔开、实为同一 NB-ARC 域；而 RPP2A 的 `TIR(3) NBARC(5) TIR(2) NBARC(5) LRR`
  是真串联。仅靠位置无法区分，改用「块规模 + 后随 NB-ARC」双条件后两者都能正确折叠。
- 验证（RefPlantNLR 400 条）：
  - DNA 默认档架构一致率 **356/399 → 361/399（90.48%）**，修好 5 例（CaMi、Gpa2、Mi-1、
    RPM1、Rx）、**零回退**；
  - 蛋白模式 **378/398（94.97%）保持不变**（RPP2A 的串联结构未被误伤）；
  - 新增测试：`domain_string_enforces_domain_order`（含 RPM1 形态、RPP2A 形态、单 motif CC
    仍有效、串联 NB-ARC、无 NB-ARC 五种情形）；`nlr-domain` 中「CC 位于 NB-ARC 之后」的断言
    由 `NBARC-CC` 更新为 `NBARC`（行为变更，已在测试内注明理由）；
  - 全量测试 67 项通过。

### 结构域判别优化 Phase 2：DNA 模式融合 HMM 证据（2026-10-12）

- 目标：A 类误差（漏判 N 端域，23 例）的根因是 DNA 模式**只能靠 PWM** 判 N 端域——内置
  20-motif 库没有 RPW8 型 CC，且 Rx_N 型 motif_16/17 对发散序列不敏感。实测这 23 个基因
  在 DNA 模式下 CC motif 命中数为 0/23，而蛋白 HMM 在 20/23 上检出了 N 端域。
- 关键设计数据（276 个真值含 CC 的基因）：PWM 灵敏度 89.9%、HMM 94.6%、
  **并集 97.8%**、交集 86.6%；非 CC 基因上 PWM 误报 9 例、HMM 误报 0 例。
  结论：N 端域必须取**并集**，单靠任一检测器都不够——这也解释了为何 DNA 与蛋白两种模式
  对同一基因的判定有 11.8% 不一致。
- 改动：
  - `nlr-core`：`MotifList` 新增 `hmm_nterm: Option<DomainCategory>` 字段（`new()` 默认
    `None`，不破坏既有构造点），`domain_string()` 在有值时将其置于域串最前。
  - `nlr-cli`：新增 `annotate_hmm_nterm()`——对每个位点取 `[span.start-1500, span.end+300]`
    窗口，六框翻译后跑内置 HMM，取 NB-ARC bitscore 最高的读码框，用与蛋白模式**同一个**
    `nlr_domain::classify()` 判定，取其 N 端域。**只在 PWM 完全没有 CC/TIR 命中时**写入，
    保证只做加法。
  - `nlr-cli/src/main.rs`：`write_outputs` 改为 `&mut RunResult`，基因组只读一次
    （同时供 HMM 补域与位点序列提取，避免重复读取），并在写任何输出前完成补域。
- 验证（RefPlantNLR 400 条）：
  - DNA 默认档架构一致率 **361/399 → 382/399（95.74%）**，修好 21 例、零回退；
  - 蛋白模式输出**逐字节不变**（该路径不经过 `write_outputs`）；
  - 耗时：小数据集 0.1 s → 0.89 s；TAIR12 真实基因组 **39.32 s，与 Phase 2 前持平**
    （仅 216 个位点需要补域，开销可忽略），内存 162 → 167 MB；
  - TAIR12 位点数不变（216），类型分布向完整型偏移：CNL 30→33、NBARC-LRR 21→23、
    TIR 36→38、新增 CC-NBARC 1。
- 残余误差 17 例的构成（重要，说明瓶颈已转移）：
  - 9 例属**金标准自身标注不一致**（Structure 无 C 但 Subclass 写 CNL），见 Phase 3；
  - RPP2A 是 DNA 模式把串联结构切成 2 个位点（过度切分问题，非域判定问题）；
  - 其余为真实检测缺口：Bs2/Rpa1/Ptr1（HMM 与 PWM 均未检出 CC）、Prv^1/Prv^2
    （尾随 NB-ARC 未检出）、Pb1（PWM 误出 TIR 致 N 端判成 TIR）、YrSP（漏 LRR）。

### 结构域判别优化 Phase 3：评估口径对齐（2026-10-12）

- 问题：原评估只看「整串完全匹配」，把「少一个域」与「多一个域」同等计错，且被金标准
  自身的矛盾严重干扰。`Structure`（域字母，来自 Pfam/Gene3D/SSF 命中）与 `Subclass`
  在部分条目上互相矛盾——例如 CD3.1 的 Subclass 写 `NLR;CNL`、Structure 却是 `NL`。
- 改动（`benchmark/Script5-Run-GoldNLR-Accuracy-Eval-V1.R` 新增第 5 节）：
  - 逐域指标 `accuracy_by_domain.csv`：CC / TIR / NB-ARC / LRR 各自算 TP/FP/FN/TN 与
    灵敏度、精确率、特异度、F1（CC/TIR 只认 N 端位置，NB-ARC/LRR 只要出现）；
  - 歧义标记：`Subclass` 含 CNL 但 `Structure` 无 C/B 的条目标为 `ambiguous`，
    主指标同时给出「含歧义」与「仅无歧义」两套数字。
- 结果（Phase 1+2 之后，RefPlantNLR 400 条）：
  - 召回：dna_default / dna_flexible / dna_extended 均 399/400（99.75%，仅漏 Wed），
    protein 398/400（99.50%，漏 CW9、Pb1）；
  - 架构一致率：dna_default **95.74%**（改造前 89.22%）、dna_extended 88.22%
    （改造前 82.66%）、protein 94.97%；
  - **剔除 12 条歧义条目后 dna_default 达 97.42%**；
  - 逐域 F1（dna_default）：CC 0.989、TIR 0.995、NB-ARC 1.000、LRR 0.999；
  - 12 条歧义条目：CD3.1、RPG5、NtTPN1、Pi54、Rpi-amr3i、Lr1、Lr21、Tsn1、Pm60、
    Pm60a、Pm60b、PmR1（17 例残余误差中有 7 例落在此桶）。
- 未采纳项：原计划在 `.nlr.txt` 增加 `type`（TNL/CNL/RNL/NL/TN）列——会破坏既有 7 列
  格式并影响 Script5/Script7 等下游解析，收益有限，暂缓（层级类型改在评估脚本侧派生）。
- 归档同步：`benchmark/arabidopsis_tair12/outputs/fastnlr_dna_default/` 用当前软件重跑
  （216 位点 / 145 完整 / 39.4 s / 167 MB），`Fig11–Fig13` 与模块 README 同步更新。

### 组装层攻坚：吸收过度切分的冗余位点（2026-10-12）

- 问题：金标准 402 个位点里有 3 个基因被切分（Prv^1、Prv^2、RPP2A），TAIR12 上有 2 对
  自身重叠/嵌套位点。根因是 `merge_seeds` 的 **rank 单调约束**与 `elongate` 的**只允许 rank
  递增**都无法处理**串联重复结构**：RPP2A 的 `TIR-NBARC-TIR-NBARC-LRR` 第二单元 TIR 的
  rank(1~3) 低于前一段 NB-ARC(4~12)；Prv 的 `TIR-NBARC-LRR-NBARC` 尾随 NB-ARC rank(4)
  低于 LRR(14)。残留 motif 各自组装，形成嵌套位点或下游孤立小位点。
- 改动（`crates/nlr-assemble/src/lib.rs` 新增 `absorb_redundant_loci()`，在 `elongate` 之后、
  命名之前调用，与结构域判定完全解耦）：
  1. **完全嵌套**：位点 B 跨度完全落在同链位点 A 内 → 并入 A。同链嵌套的两个独立 NLR 在
     自然界不成立，只可能是同一基因的重复单元未合并；
  2. **下游孤立 NB-ARC 段**：B 只含 NB-ARC/linker/NA 类 motif，且落在某个「以 LRR 结尾」的
     位点 A 下游 `distance_for_elongating`(2500 bp) 之内 → 并入 A（串联的尾随 NB-ARC）。
  两类都只并入 motif 并丢弃 B，不改变 A 的边界。
- 实现前后的数据验证（先 Python 原型再写 Rust）：金标准 402 → 399 位点，
  **过度切分 3 → 0**，架构一致率 +3 例且零回退；TAIR12 216 → 214 位点，
  **自身重叠/嵌套 2 → 0**；运行时间与内存均无变化（39.4 s / 167 MB）。
- 回归测试 2 项（`absorb_nested_tandem_repeat_unit`、`absorb_trailing_nbarc_segment`），
  已实测在停用吸收时**必然失败**，证明测试有效。

### 输出精简：移除 motif 编号清单（2026-10-12）

- 需求：输出结果文件中不再显示 motif 编号清单。
- 改动：
  - `crates/nlr-output/src/lib.rs`：`.nlr.txt` 由 7 列改为 **6 列**
    （seqname / name / domain-class / start / end / strand）；`.loci.fasta` 的头部描述去掉
    motif 清单，保留 `>名称 基因组 起-止 strand:±`。
  - `crates/nlr-core/src/motif_list.rs`：删除已无调用方的 `motif_list_string()`。
  - 未受影响：`.nlr.gff`（属性只有 `name=…;nlrClass=…`）、`.nlr.bed`、`.motifs.bed`
    （命中明细与坐标完整保留，因此信息并未丢失，只是不再冗余到 txt）。
  - 依赖方同步：`Script5` / `Script7` 的列名改为 6 列；`Script7` 的完整度判定改用
    `.nbarc.fasta`（含 P-loop 的位点）+ 类别串含 "LRR"，精确复现 `is_complete_nlr`；
    两处测试断言同步更新。
- 顺带查实一个**既有口径不一致**：`.nlr.bed` 的 itemRgb 不能用来判完整度——"含终止子"红色
  优先级高于完整绿色，TAIR12 上绿色 143 而摘要 complete 为 145（差的 2 个是完整但含终止子）。
  已在 Script7 注释中记录该陷阱。
- 最终评估（RefPlantNLR 400 条，逐域指标见 Phase 3 一节）：
  - 召回：dna_default / flexible / extended 均 399/400（99.75%），protein 398/400（99.50%）；
  - 架构一致率：dna_default **96.49%**（最初 89.22%）、dna_extended 95.24%、protein 94.97%；
  - **剔除 12 条金标准歧义条目后 dna_default 达 98.19%**；
  - 逐域 F1（dna_default）：CC 0.989、TIR 0.995、NB-ARC 1.000、LRR 0.999；
  - 过度切分 0 个基因。
- 归档同步：TAIR12 模块输出重跑（214 位点 / 145 完整 / 39.4 s / 167 MB），
  `Fig11–Fig13`、模块 README 同步更新。

### 文档精简、v1.5.0 发版与一次 CI 事故的处理（2026-10-12）

- 文档重构：README 由 386 行精简到 174 行（保留简介、安装、快速开始、特性概要、性能与链接），
  把完整参数表、输出格式、扩展 motif 库等细节迁移到官网主页 `docs/index.html`；主页新增
  "CLI reference" 与 "Output formats" 两个使用说明模块（原 Usage 区块扩展为使用说明），
  新增 "Extended RNL / TIR library" 小节，特性卡片补齐柔性播种/HMM 域判定/串联吸收三项，
  并调整各区块底色 class 维持浅深交替。版本号全量升到 **1.5.0**（含 Cargo.toml/Cargo.lock、
  各 crate 依赖约束、conda 配方、README 徽章与 GFF 示例、官网徽标与下载链接）。
- **CI 事故与修复**：首次推送后 CI 失败，原因是新增的日志时间戳测试在 CI 上解析错位——
  该环境的 `NO_COLOR` 未设置，tracing-subscriber 默认着色，日志行首带 ANSI 转义码
  （本机因 `NO_COLOR=1` 而未暴露）。这同时暴露一个真实缺陷：**stdout 被重定向到文件或管道时
  日志仍会写入 ANSI 转义码**，会污染日志文件并破坏下游解析。
  - 修复：`init_logging` 改为 `with_ansi(stdout().is_terminal() && NO_COLOR 未设置)`，
    仅在真实终端着色并遵循 `NO_COLOR` 约定；测试侧也增加 ANSI 清洗以提高健壮性。
  - 实测：模拟 CI 条件（`env -u NO_COLOR` + 管道输出）下日志已无转义码；
    `env -u NO_COLOR cargo test` 与常规 `cargo test` 均 69 项全通过。
  - 首次打出的 `v1.5.0` 标签指向未修复的提交（其自带测试在 CI 上会失败、且二进制仍会向管道
    写转义码）。本可重打标签，但 `gh` token 无删除 Release 的权限（403），而发布动作
    （`taiki-e/upload-rust-binary-action`）**不会覆盖**同名资产——重打标签会留下
    「标签指向新提交、Release 资产却是旧构建」的不一致状态。因此改为**增量发补丁版 v1.5.1**
    （历史只追加、无需网页端操作）：版本号在 Cargo.toml/Cargo.lock、各 crate 依赖约束、
    conda 配方、README 徽章与下载链接、官网徽标与 GFF 示例处统一升到 1.5.1。
    `v1.5.0` 的 Release 仍保留（二进制可用，仅缺日志着色修复），如不需要可在网页端删除。

## 2026-10-14

### benchmark 目录重整：功能分层 + 流水线编号（零删除）

- 动机：`benchmark/` 经过性能测评（v1.4.0 → v1.5.x 两代）与两条准确性流水线后，出现多代产物
  混杂——`outputs/`（v1.4.0）与 `outputs_v2/`（120 次运行 / 2.0 GB）并存、`metrics/` 里
  `goldnlr_*` 与 `accuracy_*` 重复、`figures/` 里 Fig1–10 是 FigA*/FigB* 的旧代副本、
  `Script3/4` 已标 SUPERSEDED 仍与现行脚本平铺，另有一处误建的嵌套空目录 `benchmark/benchmark/`。
- 新结构（`benchmark/README.md` 为唯一协议文档）：
  `scripts/`（S1–S8 按流水线编号）、`data/{RefPlantNLR,TAIR12}`、`outputs/{performance,
  refplantnlr_accuracy,tair12}`、`metrics/`、`figures/`、`logs/`、`archive/`。
- 脚本重命名（AGENTS.md 的「序号-功能-版本.编程语言」风格）：
  `Script1→S1-Run-Performance-Benchmark-V1.R`、`Script2→S2-Plot-Performance-Benchmark-V1.R`、
  `Script5→S3-Run-RefPlantNLR-Accuracy-V1.R`、`Script6→S4-Plot-RefPlantNLR-Accuracy-V1.R`、
  `Script7→S8-Plot-TAIR12-Locus-Profile-V1.R`；新增 `S5`（建金标准）、`S6`（运行+评测）、
  `S7`（对比图）。`Script3/4` 移入 `archive/scripts/`。
- 指标去重命名：`accuracy_*` → `refplantnlr_{gene_table,summary,by_class,by_type,by_structure,
  by_domain}.csv`；性能三个文件名保持不变（S1 的断点续跑依赖 `run_manifest.csv` /
  `genome_benchmark_raw.csv`）。S1/S2/S3/S4/S8 只改路径常量与产出文件名，运行语义不变。
- 归档（只移动不删除，清单与「可否删除」标注见 `benchmark/archive/README.md`）：
  v1.4.0 的 `outputs/`（196 文件）、旧代 Fig1–10、`metrics` 里 47 个 legacy 文件、
  `test/`（玉米冒烟）、`protein_dna_mode_misrun/`（脚本 bug 误跑残留）、
  `GoldNLR-RefPlantNLR/outputs/`（旧 CLI 接口产物）。空目录另列出（不含文件）。
- **核对**：整理前后均为 1937 个文件，逐目录计数一致，无丢失。
- 重命名后回归验证：S3 重跑（复用已有输出）指标与整理前一致（架构一致率 96.49%/98.19%、
  逐域 F1 CC 0.989 / TIR 0.995 / NBARC 1.000 / LRR 0.999）；S4 与 S2 重出图正常。
  注意 `run_manifest.csv` 中历史行仍记录旧的绝对路径（作为执行记录保留），`--rerun` /
  `--collect-only` 重新解析时需按新路径自行核对。

### 新增流水线：TAIR12 真实基因组准确性评测（S5 → S6 → S7）

- 目标：在拟南芥 TAIR12（`GCA_978657495.1`，142.5 Mb 完成图）+ NCBI 官方注释上，评估 FastNLR
  DNA 模式直接注释出的 NLR 位点是否准确——即与注释翻译蛋白经 HMM 检索得到的真实 NLR 基因位点比对。
- **金标准（S5）**：`md5 校验 → 解压 → gffread 翻译 → 断言一致性 → hmmscan --cut_ga → 蛋白汇总到基因`。
  - gffread 从「基因组 + GFF」翻译出 32,589 条蛋白，与 NCBI 官方 `protein.faa` **逐条完全一致
    （32,589/32,589）**，脚本内为硬断言，不一致即终止。
  - 用仓库内置 5 个 Pfam 模型 + 系统 hmmscan 3.3.2 独立判定（不使用 FastNLR 自身实现，
    避免循环论证）：467 条域命中 / 277 条蛋白；NB-ARC 231 蛋白 → 166 基因；TIR 180 蛋白；
    CC（Rx_N/RPW8）40 蛋白 → 27 基因；PF12061 在拟南芥无命中。
  - 金标准 = NB-ARC 基因 ∪ TIR-only 基因 = **198 个**（TNL 98 / CNL 27 / NB-only 41 / TIR-only 32；
    Chr1 66 / Chr2 8 / Chr3 18 / Chr4 38 / Chr5 68）。
  - 口径要点：基因级（任一 isoform 命中即命中）；不含 LRR 类 Pfam（Pfam LRR 在 NLR 上灵敏度差）；
    CC 灵敏度有限（27 基因）故 CNL 召回会被系统性低估，已在 README 与图注标注。
- **被评估两档（S6）**：`fastnlr_dna_default`（内置 20 motif）与 `fastnlr_dna_extended`
  （28 motif 扩展库），`-t 8 --progress off` 跑全基因组。
- **匹配判据**：主判据 = gene 跨度（GFF `gene` start–end）任意重叠 ≥1 bp；稳健性对照 =
  CDS 跨度重叠与覆盖率 ≥50%。

| 指标 | default | extended |
|---|---|---|
| NLR 位点 / 完整位点 | 214 / 145 | 252 / 145 |
| 金标准基因召回 | **194/198 = 98.0%** | 196/198 = 99.0% |
| 位点精确率 | **92.5%**（198/214）| 80.6%（203/252）|
| 未匹配位点（基因间区 / 非编码基因 / 非 NLR 蛋白编码基因）| 16（4 / 9 / 3）| 49（11 / 10 / 28）|
| 拆分基因 / 合并位点 | 12 / 10 | 21 / 18 |
| gene 跨度覆盖率中位数 | 0.68 | 0.69 |
| 类型一致性 | 86.6% | 83.7% |
| 耗时 / 峰值内存 | 39.4 s / 164 MB | 54.6 s / 162 MB |

- 分类型召回（default）：TNL 98/98、CNL 25/27、TIR-only 31/32、NB-only 40/41；
  未检出的 4 个基因为 AT4G23440(NB-only)、AT4G33300(CNL)、AT5G45090(TIR-only)、AT5G66630(CNL)。
- **结论**：默认档在 1 bp 重叠口径下召回 98.0% 的 HMM 金标准 NLR 基因、位点精确率 92.5%，
  拆分/合并各 ≤12 例，覆盖率中位数 0.68；扩展库把召回提到 99.0%，但精确率降到 80.6%，
  新增位点主要压在非 NLR 蛋白编码基因上（3 → 28），且类型一致性略降、拆分/合并翻倍，
  因此**默认档是更合适的生产配置**，扩展库更适合"宁可多报不可漏报"的候选位点枚举场景。
- 与设计预期的偏差：侦察阶段预判「16 个未匹配位点里约 13 个在基因间区」，实际做细粒度归因后
  分为 4 个基因间区 + 9 个压在非编码基因上（合计 13，与预判一致，只是归因更细）。
  另外侦察阶段按蛋白级口径算出的 TIR-only 为 35，改为基因级（任一 isoform 含 NB-ARC 即归 TNL）
  后为 32——3 个基因同时拥有 NB-ARC 型与 TIR-only 型 isoform，基因级口径更严谨，金标准总数
  仍为 198（与设计一致）。
- **图件（S7）**：`Fig14-TAIR12-accuracy-genome-map.pdf`（a 双轨染色体图：上轨金标准基因
  命中/漏检、下轨 FastNLR 位点 命中/未匹配；b 各染色体 命中基因/漏检基因/未匹配位点 堆叠柱，
  两档分面）、`Fig15-TAIR12-accuracy-summary.pdf`（a 分类型召回、b 未匹配位点归因、
  c 已检出基因的覆盖率分布）。描述性 Fig11–13 由 S8 用新路径重出，内容不变
  （214 位点 / 145 完整）。
- 实现中踩到的两个坑，已在脚本注释中固化：
  1. **gffread 0.12.7 不能读 `.gz`**：基因组与 GFF 必须先解压成 `raw/*.fna` / `raw/*.gff`
     （该约定原本就写在数据集 README 里），否则报 "sequence lines in a FASTA record must have
     the same length"（FASTA）或 "invalid start coordinate"（GFF，随后栈溢出）。
  2. **`--seed-combination` / `--signature` 的值含分号**：`system2` 经 shell 传递，不 `shQuote`
     会被拆成多条命令，实测得到一份"参数被截断"的假结果（186 位点），因此 S6 对三个长串参数
     统一 `shQuote`，并在 README 的运行手册中保留该写法。
  3. **data.table 分组汇总里不能用分组列做下标**：`by = .(config, gold_class)` 内 `gold_class`
     退化为单值，`gold_class[detected]` 会静默产出全 NA，导致 `tair12_by_type.csv` 的
     `type_agreement` 整列为 NA。改为先在 `gene_tab` 上算出独立的 `type_match` 列再求均值。
     修复后分类型类型一致性（default）：TNL 0.98、CNL 0.92、TIR-only 0.81、NB-only 0.60
     （NB-only 偏低是因为部分位点把邻近域的 TIR/CC 一并纳入，导致位点类别首域不是 NBARC）。
- 文档：新增 `benchmark/README.md`（目录地图 / 数据集 / 三条流水线 / 口径定义 / 运行手册 /
  指标文件 / 结果摘要 / 图件索引 / 已知事项）、`benchmark/archive/README.md`（归档清单与
  文件计数核对）、`benchmark/data/RefPlantNLR/README.md`；重写 `benchmark/data/TAIR12/README.md`
  （新路径、金标准口径与结果、两档基线与结论、环境核查）。
- 说明：`benchmark/` 仍在 `.gitignore` 内，本次整理与产出一律只留本地；`.gitignore` 未改动。

### TAIR12 误差诊断：漏检基因与未匹配位点的逐条归因（2026-10-14）

- 动机：S6 的评测显示 default 档漏检 4 个金标准基因、有 16 个位点不与任何金标准基因重叠。
  这两类误差需要用独立证据说清楚"到底是什么原因"，否则召回率/假阳性数字会被误读。
- 新增脚本 `benchmark/scripts/S9-Diagnose-TAIR12-Errors-V1.R`（只读 S5/S6 产出，不改现有表）：
  - **Miss 三层证据**：全蛋白组 `hmmscan -E 10`（无 GA）复核打分 → 基因内/±5 kb 的 accepted
    motif 命中 → 基因 ±20 kb 局部重跑四档（default / relaxed / sensitive / extended）。
  - **Unmatched 三路复核**：重叠基因 CDS 蛋白 / 位点六框 ORF / P-loop 锚定片段，合并成一个
    FASTA 后一次跑 `--cut_ga` 与 `-E 1e-3`，判定 `gold_standard_miss` / `divergent_partial` /
    `fragment_call` / `unresolved_spurious`。
  - **口径敏感性**：S5 `--relaxed` 新增输出 `tair12_gold_standard_relaxed.csv`（宽口径
    E≤1e-3，317 个基因，含坐标），S9 用它重算召回区间。
- 结论：
  - **漏检 4 个基因全是真实漏检**：AT4G23440(TN17)、AT4G33300(ADR1-L1)、AT5G66630(DAR5)
    基因内 accepted motif 命中为 0（卡 `p<1e-5` 扫描阈值，放宽到 1e-3 即可召回）；
    AT5G45090 基因内有 1 条 motif_18（p=2.5e-9）但单 TIR 不构成播种组合（卡播种）。
    局部重跑：sensitive 4/4 召回、extended 2/4、relaxed(--no-flexible-seed --relaxed-seed 3) 0/4。
  - **16 个未匹配位点里只有 11 个是工具自身产物**：5 个含 P-loop 的大位点（1.9–4.0 kb）经六框
    ORF 复核在 GA 阈值下命中 NB-ARC（最高 200 bits）或 TIR（最高 179 bits），重叠的注释基因
    在 TAIR12 里是 other/lncRNA 或被标成 "putative ATP hydrolase"，属**金标准漏判**；
    11 个是 151–619 bp、无 P-loop 的片段型命中（5 个含终止子）。extended 档 49 个未匹配里
    45 个是片段型。
  - **有效指标**：default 精确率 92.5% → 有效 94.7%（下界）/94.9%（上界），有效假阳性 11 个
    （5.1%）；extended 80.6% → 81.9%/82.1%。召回 98.0% 无需修正。宽口径金标准（317 基因）下
    召回降到 61.2%/63.1%，下降源于金标准分母膨胀（两档检出的绝对基因数仅 194/200）。
- 实现中修掉的三个坑：
  1. `--relaxed-seed` 在 v1.5.x 与默认开启的柔性播种**互斥**（直接报错），relaxed 档必须写成
     `--no-flexible-seed --relaxed-seed 3`；
  2. `.nbarc.fasta` 是 **gap 填充的比对片段**（含 `-`），HMMER 拒绝解析，需先剔除 gap；
  3. hmmscan 失败时会**留下残缺的 domtbl**，导致后续运行误用缓存——改为写 `.tmp` 成功后改名。
  另外 `fread` 读 BED 必须显式 `sep="\t"`，否则自动分隔符会把 12 列解析成 65 列。
- 脚本内置自校验（漏检归因必须复现 3×scan_threshold + 1×assembly_seed、nlr49/nlr14 必须判为
  gold_standard_miss、sensitive 必须召回全部 4 个漏检基因），不满足即报错退出。
- 产物：`benchmark/reports/TAIR12-error-diagnosis-V1.md`（完整诊断报告）、
  `metrics/tair12_diagnosis_{miss,locus,stage,summary}.csv`、`metrics/tair12_gold_standard_relaxed.csv`、
  `figures/Fig16-TAIR12-error-attribution.pdf`、`outputs/tair12/diagnosis/`；
  `benchmark/README.md` 新增「八、误差诊断」并把有效指标并入「六、结果摘要」。

## 2026-10-15

### 项目体检、两处文档/CLI 瑕疵修复与 v1.5.2 发版

- **全量体检结论（健康）**：`cargo test --workspace --locked` 69/69 通过、`cargo fmt --check` 与
  `cargo clippy -D warnings` 无告警、`--locked` 通过（依赖锁一致）、二进制与源码同步；代码中无
  `TODO`/`FIXME`/`unimplemented!` 残留。文档侧：版本号在 Cargo.toml / conda 配方 / README 徽章与
  下载链接 / 官网四处一致；`fastnlr --help` 的参数集与官网 CLI reference 逐条对应；现行文档中
  已无旧脚本名与旧路径残留。数据侧：`genome/` 18 物种 14 GB，TAIR12 的 MD5 每轮由 S5 自动校验。
  结果侧：README 引用的性能四组数字、RefPlantNLR 指标、TAIR12 指标与误差诊断数字均与
  `benchmark/metrics/*.csv` 逐条核对一致。
- **修复两处小瑕疵**：
  1. `--relaxed-seed` 与默认开启的柔性播种互斥这件事此前只写在官网表格里，`fastnlr --help`
     未提示，用户会直接撞运行时报错——已在参数帮助文本中补上
     "Mutually exclusive with the default flexible seeding — pass --no-flexible-seed together with it"。
  2. 文档中 TAIR12 的耗时/内存是首跑的 39.4 s/164 MB 与 54.6 s/162 MB，与复跑后的
     `summary.txt`（39.9 s/164 MB、54.7 s/158 MB）有毫秒级漂移——已按当前产物对齐
     （`benchmark/README.md`、`benchmark/data/TAIR12/README.md`、诊断报告）。
- **v1.5.2 发版**：版本号在 Cargo.toml / Cargo.lock（10 个工作区包）、各 crate 依赖约束、
  README 徽章与下载链接、官网（徽标、二进制下载地址、GFF 示例的 `##source-version`）处统一升到
  1.5.2；推送 `main` 后打 annotated tag `v1.5.2`，Release 工作流三平台构建全部成功
  （Linux / Windows / macOS 二进制 + sha256），CI 与 GitHub Pages 部署同时通过。
- **发布校验**：从 Releases 下载 Linux 资产，sha256 与官方 `.sha256` 完全一致（`422230c1…`）；
  解包后二进制自报 `fastnlr 1.5.2`，且包含本次新增的 `--relaxed-seed` 互斥提示。
- **conda 配方同步**：source sha256 更新为 v1.5.2 标签归档的校验和
  `b438f737e5d6d840710160647bb3cc944cd63476044b5a6877a9843ce0b2964b`（连续两次下载校验一致），
  与 v1.5.0 / v1.5.1 的发版流程保持一致。

### FastNLR vs NLR-Annotator v2.1b 横向评测（TAIR12，128 线程）

- 新增 `benchmark/scripts/S10-Run-Tool-Comparison-V1.R`（运行 + 评测）与
  `S11-Plot-Tool-Comparison-V1.R`（Fig17 / Fig18），报告
  `benchmark/reports/TAIR12-tool-comparison-V1.md`。
- **对照设计已证实为最干净形态**：`/home/zjw/software/mot.txt` 与 `store.txt` 和 FastNLR
  **内置库 md5 完全一致**（`c307517d…` / `0d1ede50…`，20 motif），`-t/-n` 语义与输出格式同源，
  因此差异只来自实现与流水线阶段。两工具统一读**未压缩** fasta（NLR-Annotator 不支持 gz）、
  同为 `-t 128`、`-n` 默认 1000，各 3 次重复、串行执行、每次前预热页缓存。
- **速度与内存**（TAIR12 142.5 Mb，128 线程，3 次中位数）：FastNLR 默认档 **6.01 s / 164 MB**
  （5.27/6.09/6.01 s），扩展库档 7.96 s / 165 MB；NLR-Annotator **465.2 s / 2 266 MB**
  （465.2/458.8/526.6 s，第 3 次峰值 4 316 MB）。→ **速度 77×、内存 14×**。
  CPU 利用率 FastNLR 约 86×，NLR-Annotator 仅 **9.4×**——同参数下其 JVM 线程池未有效扩展。
- **准确性**（198 个 HMM 金标准基因）：FastNLR 召回 **97.98%**（194/198）、假阴性率 2.02%、
  精确率 92.52%、假阳性率 7.48%；NLR-Annotator 召回 **92.93%**（184/198）、假阴性率 7.07%、
  精确率 94.42%、假阳性率 5.58%。按 S9 的证据修正后有效精确率 94.74% vs 96.88%、
  有效假阳性率 5.14% vs 3.05%。
- **召回集是严格超集**：NLR-Annotator 检出的 184 个基因 FastNLR 全部检出，且多检出 10 个
  （ADR1、ADR1-L2、ADR1-L3、NRG1.1 四个 RNL/helper + 5 个 TIR-only 与 1 个 NB-only），
  正是 S9 归因的
  "柔性播种 + HMM 补 N 端域"生效的场景；两者共同漏掉的仍是那 4 个（TN17、ADR1-L1、DAR5、
  AT5G45090），瓶颈在 20-motif PWM 库的检测下限。
- **结果一致性**：NLR-Annotator 的 197 个位点 **100%** 与 FastNLR 位点重叠，其中 **186 个
  （94.4%）坐标完全相同**，Jaccard 0.921；FastNLR 独有 17 个位点（TIR 11、CC-NBARC-LRR 3、
  NBARC-LRR 2、NBARC 1，其中 12 个命中金标准基因、5 个为片段），NLR-Annotator 独有 0 个。
- 实现中修掉的两个坑：`data.table(key = ...)` 会把 `key` 当保留参数（列名改 `metric`）；
  `chr2acc[match(seqname, ...)]` 方向写反导致染色体映射错（应为 `chr_tab$chr[match(...)]`）。
  另发现 nohup 后台任务会随 exec 会话结束被杀，改用 `setsid` 才能跑完 23 分钟的 Java 任务。
- 图件：`Fig17-TAIR12-tool-comparison.pdf`（速度 / 内存 / 准确性指标 / 位点集合四面板）、
  `Fig18-TAIR12-tool-three-track.pdf`（金标准 / FastNLR / NLR-Annotator 三轨染色体图）。
  `benchmark/README.md` 新增「九、工具横向对比」并把结论并入「六、结果摘要」。

### 修复：柔性播种接受"同一 motif 重复"造成的退化种子

- **问题定位**：与 NLR-Annotator 的横向对比显示 FastNLR 召回率明显更高（97.98% vs 92.93%），
  但位点级假阳性率略高（7.48% vs 5.58%）。逐条归因后，16 个未匹配位点里 5 个是**金标准漏判**
  （两工具都有），10 个是片段型命中，其中 3 个来自刻意加入的 `[18,15]`/`[15,13]` 双 TIR 种子
  （这套种子同时贡献了 5 个 TIR-only 基因的召回，是净收益为正的既有取舍，不动），
  **只有 1 个是真正的缺陷**：`OZ408684.1_nlr2`（Chr2:11,581,849–11,582,467，基因间区，
  619 bp、无 P-loop、无 LRR）。
- **根因**：柔性播种判据原为"NB-ARC motif 数 ≥ 3 且 rank **单调不减**"。8 个 NB-ARC motif 的
  rank 两两不同（4/5/6/7/8/9/10/11），因此"不减"与"严格递增"只在**同一 motif 被重复命中**
  时才有差别——该位点上连续命中三次 `motif_2`（rank 11,11,11），被当成"按 rank 排列"而播种。
  这类重复不构成结构域的先后顺序。
- **修复**（`crates/nlr-core/src/signature_def.rs`）：
  - `is_seed_flexible` 判据由 `rank(a) <= rank(b)` 收紧为 `rank(a) < rank(b)`（严格递增）；
  - **片段预过滤保持宽松判据**（新增 `is_seed_flexible_loose`，`has_signature_flexible` 用它）。
    这一步是关键：预过滤面对的是**整个片段的命中串**，若一并收紧，片段尾部的一个重复 motif
    会让该片段里合法的种子前缀被整段丢弃——首次尝试时正是这样误删了
    `OZ408687.1 1,365,557–1,367,650` 的 CC-NBARC-LRR 位点（motif 1,4,10,1，丢掉 1 个基因召回）。
  - 回归测试：3 次同 motif、以及"合法前缀 + 尾部重复"都必须**不能**播种，
    但 `is_seed_flexible_loose`/`has_signature_flexible` 必须继续放行。
- **验证**（全部用修复后的 release 二进制重跑）：
  - 单元/集成测试 69/69 通过，`cargo fmt --check` 与 `cargo clippy -D warnings` 干净；
  - TAIR12：位点 **214 → 213**（只少了那一个基因间区位点），召回 **不变 194/198 = 97.98%**，
    精确率 **92.52% → 92.96%**，位点级假阳性率 **7.48% → 7.04%**，有效假阳性率 5.14% → 4.69%；
  - RefPlantNLR 四档金标准：召回 399/400（protein 398/400）、架构一致率 96.49%/95.24%/94.97%
    **逐项不变**；
  - 多个物种位点数：拟南芥 halleri 332、lyrata 275 不变；黄瓜 77（其 `--no-flexible-seed`
    为 74，与 v1.4.0 基线一致，说明多出的 3 个是 v1.5 柔性播种的既有行为，非本次修复引入）。
- **文档同步**：`crates/nlr-cli/src/main.rs` 的 `--flexible-seed` 帮助文本、`crates/nlr-core/README.md`、
  `crates/nlr-assemble/README.md`、`docs/index.html` 的 CLI 表统一改为"strictly increasing rank"；
  `docs/flexible-seeding-upgrade.md` 增加 2026-10-15 修订说明（含预过滤为何保持宽松）。
- **benchmark 侧刷新**：S6/S7/S8/S9/S10/S11 全部用修复后的二进制重跑，
  `metrics/*`、`Fig11–Fig18` 与两份报告的数字同步更新；
  另修掉 S10 自身的一个键名 bug（`stops_by_cfg` 用 `names(CFG)` 命名但取值用 `cfg$name`，
  导致 `has_stop` 恒为 FALSE）。

## 2026-10-16

### v1.5.3 发版：柔性播种退化种子的修复

- 版本号在 Cargo.toml / Cargo.lock（10 个工作区包）、各 crate 依赖约束、README 徽章与下载链接、
  官网（徽标、二进制下载地址、GFF 示例的 `##source-version`）统一升到 **1.5.3**，
  conda 配方的 `version` 同步（source sha256 在打 tag 后单独提交更新）。
- 本次发布的唯一行为变更：柔性播种判据由 rank **单调不减** 收紧为 **严格递增**，
  关掉"同一个 NB-ARC motif 连续重复"形成的退化种子（TAIR12 上基因间区的
  `OZ408684.1_nlr2`）；片段预过滤保持宽松判据不变（原因见 2026-10-15 那条记录）。
  实测影响：TAIR12 位点 214 → 213、召回不变（194/198）、位点级假阳性率 7.48% → 7.04%；
  RefPlantNLR 四档金标准指标逐项不变。
- 发布前质量门：`cargo test --workspace --locked` 69/69 通过、`cargo fmt --check` 干净、
  `cargo clippy --workspace --all-targets -D warnings` 无告警、`--locked` 构建一致。
- 发布流程与 v1.5.2 一致：推送 `main` → 打 annotated tag `v1.5.3` → Release 工作流三平台构建
  （Linux/Windows/macOS 二进制 + sha256）→ 下载 Linux 资产核对 sha256 与 `--version` →
  用标签归档的 sha256 更新 conda 配方。CI 与 GitHub Pages 同步部署。
- **发布校验**：Release（3 分 11 秒）、CI、GitHub Pages 三个工作流全部 success；
  从 Releases 下载 Linux 资产，sha256 与官方 `.sha256` 完全一致（`7da8049c…`），
  解包后二进制自报 `fastnlr 1.5.3`，且 `--flexible-seed` 的帮助文本已是
  "strictly increasing rank order (a repeated motif does not qualify)"。
- **conda 配方同步**：source sha256 更新为 v1.5.3 标签归档的校验和
  `11d6b6537abedf69afa548076c7d97cd48097ac356d6623f4ebad1c8b1d40aef`（连续两次下载校验一致）。

## 2026-10-19

### TAIR12 差异位点结构分析与可视化

- **新增脚本**：
  `benchmark/scripts/S12-Build-TAIR12-Differential-Locus-Structure-V1.R`（差异位点目录、GFF 基因模型、
  HMM 域坐标、六框 ORF/motif 轨道）与
  `benchmark/scripts/S13-Plot-TAIR12-Differential-Locus-Structure-V1.R`（主图、补充材料、中文报告）。
  两个脚本均从仓库根目录运行，复用 S5/S6/S9 与工具对比输出，未重跑全基因组扫描，也未修改 FastNLR 核心代码。
- **案例口径**：FastNLR default 漏检 4 个金标准基因；NLR-Annotator 漏检 14 个（4 个共同漏检 +
  10 个 NLR-Annotator 独有漏检）；FastNLR default 未命中金标准 NLR 基因的位点 15 个（5 个
  `gold_standard_miss` + 10 个 `fragment_call`）。补充材料共 29 个唯一案例，4 个共同漏检只在
  Section A 展示，Section B 只放 10 个独有漏检，避免重复页面。
- **结构分析结果**：
  - FastNLR 的 4 个漏检中，3 个卡在 motif 扫描阈值（AT4G23440、AT4G33300、AT5G66630），
    1 个卡在单 TIR 无法播种（AT5G45090）；extended 库恢复 AT4G33300 与 AT5G66630。
  - NLR-Annotator 的 10 个独有漏检集中在 RNL/helper（ADR1、ADR1-L2、NRG1.1）、TIR-only 与
    NB-only；这些基因全部能被 FastNLR default 检出，结构上对应柔性播种和 HMM 补 N 端域。
  - FastNLR 未匹配的 15 个位点中，5 个含 P-loop 且 ORF HMM 在 GA 阈值下命中 NB-ARC/TIR，
    形态上是完整 NLR（`gold_standard_miss`）；10 个为 151–570 bp、无 P-loop 的片段型命中
    （5 个含终止子），注释上下文为基因间区 3、非编码基因 9、非 NLR 蛋白编码基因 3。
- **图件与证据表**：
  - `benchmark/figures/Fig19-TAIR12-differential-locus-structure.pdf`：5 个代表案例，每个案例
    四层（genome/annotation、six-frame translation、FastNLR、NLR-Annotator），1 页。
  - `benchmark/figures/FigS19-TAIR12-differential-locus-structure-all.pdf`：29 页补充材料。
  - `benchmark/metrics/tair12_diff_{cases,gene_catalog,locus_catalog,gene_models,tool_loci,hmm_domains,sixframe_tracks,summary}.csv`
    与 `benchmark/outputs/tair12/differential_loci/regions.fa`。
  - `benchmark/reports/TAIR12-differential-locus-structure-V1.md`：方法口径、逐案例结构解读、
    结构特征汇总与局限性。
- **验证**：S12 内置断言 4/14/4/10、15/5/10、3/3/9、29 个案例；`nlr49` 必须含 P-loop 且有
  GA 级 NB-ARC/TIR 证据，`nlr64` 必须为 169 bp 无 P-loop 片段，`AT4G33300` 必须被 extended
  恢复，`AT5G45090` 必须归因为 `assembly_seed`。S13 输出主图 1 页、补充材料 29 页，PDF 字体为
  嵌入的 Liberation Sans，图内无中文。`benchmark/README.md` 已同步 Part 5、S12/S13 运行命令、
  指标表与 Fig19/FigS19 索引。

## 2026-10-21

### TAIR12 注释验证用 10 个 RNA-seq 样本下载

- **下载目录**：`benchmark/data/TAIR12/rnaseq/candidate10/`。包含 10 个 Col-0 样本、
  共 19 个 FASTQ 文件（9 个双端样本 + 1 个单端样本）。
- **样本覆盖**：幼苗子叶、成熟根、莲座叶第 7 片、茎生叶、茎节间、萌发后 10 天茎尖分生组织、
  花第 12 期、角果第 4 期、种子胚第 8 期和成熟花粉。
- **数据来源**：9 个样本来自 E-MTAB-7978 / PRJEB32665（Mergner et al. 2020 组织图谱），
  茎尖时间序列样本来自 E-MTAB-4202 / PRJNA268115。
- **数据量**：压缩 FASTQ 总计 **29.70 GB**。此前按 SRA run size 估算的 16.8 GB 低估了双端
  R1/R2 总量，已按 FASTQ `Content-Length` 修正。
- **下载与校验**：全部 19 个文件的文件大小与 ENA 官方 `fastq_bytes` 一致，并通过 ENA/NCBI
  官方 FASTQ MD5 校验。初次下载发现 `ERR3333396_1`、`ERR3333419_1`、`ERR3333390_1` 存在
  CRC/长度错误；旧文件已移动到 `checksums/quarantine/` 保留，随后按官方 MD5 重新下载并通过校验。
- **清单与日志**：
  - `candidate10_manifest.csv`：样本、run、组织/时期、FASTQ 路径与总大小。
  - `fastq_urls.txt` / `expected_sizes.tsv`：下载 URL 与预期大小。
  - `checksums/official_fastq.md5` / `checksums/fastq.md5`：官方与本地 MD5。
  - `checksums/final_md5_check.log` / `checksums/final_summary.tsv`：最终校验结果。
  - `logs/`：aria2c 下载与重试日志。
- **当前状态**：RNA-seq 原始数据已可用于 TAIR12 基因组回帖、转录本组装和基因模型/假基因复核；
  尚未执行比对和注释评估。

### TAIR12 RNA-seq 数据质控评估

- **检查对象**：`benchmark/data/TAIR12/rnaseq/candidate10/fastq/` 下 10 个样本、19 个 FASTQ。
- **完整性**：全部文件通过 ENA/NCBI 官方 MD5；`seqkit stats` 确认 9 个双端样本 R1/R2
  read 数逐对相同，8 个样本为 2×76 bp，1 个为 50 bp 单端。
- **质量评估**：使用 `work` 环境的 `fastp v1.0.1`，每个样本抽取前 100 万 reads/pairs，
  关闭修剪与过滤评估原始质量；另跑一次仅接头去除评估。
- **结果**：Q20 95.40%–99.00%，Q30 92.60%–96.85%，末端平均质量 Q33.3–Q38.2；
  GC 43.81%–50.32%；接头 reads 0.22%–0.87%；双端插入片段峰值 116–121 bp。
- **重复率**：多数样本 0.82%–8.48%；胚第 8 期 41.28%、茎尖单端 27.33%，更可能反映
  高表达组织转录本及 PCR 重复，不支持在 FASTQ 层全局去重。
- **判断**：数据质量足以直接进入比对/组装。无需强制质量剪切、polyG 修剪或全局去重；
  建议仅做一次轻量接头去除和质量/长度过滤，或直接使用原始 reads 比对。正式比对后应再
  检查唯一比对率、rRNA/细胞器 reads、覆盖度和剪接一致性。
- **产物**：`candidate10/qc/QC_summary.md`、`qc/fastp_qc_summary.tsv`、
  `qc/seqkit_stats.tsv`、`qc/fastp/` 与 `qc/fastp_adapter/`。

### TAIR12 RNA-seq 冒烟比对与 Fig20

- **样本选择**：使用 E-MTAB-7978 的 cauline leaf 样本 `ERR3333391`（45,325,673 对 2×76 bp reads）。
  该样本在 E-MTAB-7978 gene count 矩阵中对 `AT4G33300`、`AT5G04720` 等目标基因表达最高。
- **比对环境**：原 `work` 目录已不是有效 conda 环境，因此新建独立 `rna` 环境并安装
  HISAT2 2.2.3 + samtools 1.24，未修改 base。
- **索引与比对**：以 `GCA_978657495.1_TAIR12_genomic.fna` 构建 HISAT2 genome-only index；
  使用 `--no-unal --no-mixed --no-discordant --dta`、8 线程比对全量样本。
- **比对结果**：总比对率 **96.45%**；唯一 concordant 42,848,851 对（94.54%），多重
  concordant 869,760 对（1.92%）；BAM 为
  `benchmark/data/TAIR12/rnaseq/smoke_test/bam/ERR3333391_full.sorted.bam`
  （2.2 GB，已建索引）。
- **Fig20**：扩展 `S13` 支持可选 `--bam`，在 Annotation 与 Six-frame 之间加入 RNA-seq 轨道，
  展示窗口内 coverage、read blocks、代表性剪接弧线、最大深度和 spliced-read 数。
  使用 `ggrepel` 处理基因/结构域标签，并将主图每个工具轨道限制为一个 focal 标签，
  修复原 Fig19 中的文字重叠。
- **窗口覆盖**：5 个案例均可检出 RNA-seq 支持；例如 `AT4G33300` 窗口 34,047 条 reads、
  9,176 条剪接 reads、最大深度 1062×；`AT5G45090` 窗口 6,547 条 reads、348 条剪接 reads；
  `nlr64` 窗口仍有 538 条 reads。逐案例统计见
  `benchmark/metrics/tair12_Fig20_rnaseq_coverage.csv`。
- **产物**：
  - `benchmark/figures/Fig20-TAIR12-RNAseq-alignment.pdf`（1 页，180×340 mm）
  - `benchmark/data/TAIR12/rnaseq/smoke_test/README.md`
  - `benchmark/data/TAIR12/rnaseq/smoke_test/logs/ERR3333391_flagstat.txt`
  - `benchmark/data/TAIR12/rnaseq/reference/hisat2/TAIR12.*.ht2`

### TAIR12 多组织 RNA-seq 合并比对与 Fig20 V2

- **数据范围**：使用已通过 MD5 和质控检查的 10 个 Col-0 样本、19 个 FASTQ，
  压缩数据共 29.70 GB，覆盖根、叶、茎、子叶、茎尖、花、角果、胚和花粉。
- **比对策略**：运行 `S15-Align-TAIR12-Multitissue-RNAseq-V1.sh`，各样本以 HISAT2 2.2.3
  比对到 TAIR12 genome-only index，samtools sort/index 后合并为 pooled BAM。该方式正确
  处理 9 个双端库和 1 个单端库；逐 read 比对相互独立，因此合并 BAM 与先拼接 FASTQ
  后比对得到相同的 pooled alignment set，同时保留样本级 BAM 和日志。
- **pooled BAM**：
  `benchmark/data/TAIR12/rnaseq/alignment/pooled/TAIR12_multitissue.pooled.sorted.bam`
  约 19 GB，索引完整。flagstat 为 862,098,517 条 alignment records；primary
  672,563,251、secondary 189,535,266；properly paired 630,803,174。
- **Fig20 更新**：`S13` 的 RNA-seq 区域改为 Annotation 下方两个独立轨道：
  第一层显示 primary、MAPQ ≥10 的代表性 read exon blocks 和不同 junction 的剪接弧线；
  第二层显示完整窗口内的 IGV 风格 coverage histogram。coverage 由 samtools 在 C 层计算，
  避免 R 逐 read 展开造成的长耗时。
- **可视修正**：基因骨架增加正负链箭头；删除两工具轨道底部 2 kb 标尺；主图仅保留
  focal 基因/位点标签并将状态移至留白区；coverage 仅显示 0、中值、最大值三个 K 单位刻度，
  避免刻度文字互相遮挡。Fig20 尺寸更新为 180×390 mm，仍为单页 cairo PDF。
- **窗口结果**：`AT4G33300` 窗口 278,327 条有效比对、72,622 条剪接 read、最大深度
  7,482×；`AT5G45090` 为 43,602 / 3,238 / 1,997×；`AT5G04720` 为
  180,592 / 64,064 / 11,594×；`nlr49` 为 98,623 / 15,124 / 1,915×；
  `nlr64` 窗口仍有 31,382 条 reads 和 11,248 条剪接 read。逐案例统计见
  `benchmark/metrics/tair12_Fig20_rnaseq_coverage.csv`。
- **最终产物**：
  - `benchmark/figures/Fig20-TAIR12-RNAseq-alignment.pdf`
  - `benchmark/metrics/tair12_Fig20_rnaseq_coverage.csv`
  - `benchmark/scripts/S15-Align-TAIR12-Multitissue-RNAseq-V1.sh`
  - `benchmark/data/TAIR12/rnaseq/smoke_test/README.md`
  - `benchmark/README.md` 的 Part 6、运行命令与 Fig20 索引

### TAIR12 Fig20 单案例拆分与排版优化 V3

- **单案例输出**：Fig20 的 5 个代表案例已分别输出到
  `benchmark/figures/Fig20-TAIR12-RNAseq-alignment-panels/`，每个 PDF 独立包含
  annotation、coverage、six-frame、FastNLR、NLR-Annotator 和图例；清单见
  `Fig20_panel_manifest.csv`。
- **轨道调整**：删除 RNA-seq reads blocks/arcs 轨道，仅保留 IGV 风格 coverage；
  coverage 的 reads 与 spliced-read 统计仍写入
  `benchmark/metrics/tair12_Fig20_rnaseq_coverage.csv`。
- **注释与六框**：基因方向箭头移动到对应 3' 端；Six-frame 的 ORF 按与其
  重叠最多的 motif/HMM 类型着色，包括 NB-ARC、TIR、CC、LRR、Linker 和无 motif 灰色。
- **FastNLR 分档**：Default 与 Extended 改为上下分离的泳道组，左侧分别显示档位名称；
  Extended 轮廓改为密集点虚线，图例同步更新。
- **字体排版**：字体族参数改为 Arial，并增大标题、坐标轴、轨道标签和图例字号；
  由于系统未安装微软原版 Arial，fontconfig 实际回退并嵌入 Liberation Sans。
- **复核**：5 个单案例 PDF 均为单页、180×160 mm；总图仍为单页、180×390 mm。
  已逐图检查字号、方向箭头、coverage 刻度、六框配色、Default/Extended 标签和文字遮挡。

### TAIR12 ADR1-L1（AT4G33300）100 kb 局部漏检诊断

- **局部序列**：提取 `OZ408686.1:23,130,628-23,336,715`，长度 206,088 bp，
  覆盖 AT4G33300 基因两侧各 100 kb；同时提取局部 GFF 和代表蛋白 CAO7049742.1。
- **HMM 复核**：816 aa 蛋白在 Pfam GA 阈值下命中 NB-ARC 63.1 bits / 6.3e-21
  （aa 190-340）和 RPW8/CC 41.8 bits / 2.1e-14（aa 9-129），确认是完整
  helper NLR/RNL，而非假基因或错误注释。
- **FastNLR 五档测试**：
  - default：局部 0 个位点，基因内 0 个保留 motif。
  - `--flexible-seed 1 --require-ploop`：仍为 0 个位点，但保留 6 个真实 motif：
    `motif_1,motif_11x5`；最近 motif_1→motif_11 距离约 1.5 kb，超过 seed 内
    500 bp 检查窗口。
  - `--flexible-seed 2`：仅恢复 855 bp 的无 P-loop `NBARC` 片段，不是正确结构。
  - sensitive p=1e-3：形成 34,924 bp 的超大重复域合并位点，局部总位点数增至 36。
  - extended：命中 8 个基因内 motif，`motif_21→23→1→22→11x4`，恢复
    `CC-NBARC-LRR`；关键 seed 为相距约 402 bp 的 `motif_23→motif_1`。
- **结论**：default 漏检归因应细化为 **default motif 库/播种组合不足**，而不是
  简单的 PWM 扫描阈值不足。最优修复方向是把 RNL/helper CC motif 与 `CC,1`
  播种规则纳入 default，或增加强 P-loop 到 LRR 的结构特异性播种；不建议全局
  降到 p=1e-3。
- **产物**：
  - `benchmark/scripts/S16-Diagnose-TAIR12-ADR1L1-Recall-V1.R`
  - `benchmark/reports/TAIR12-ADR1L1-recall-diagnosis-V1.md`
  - `benchmark/outputs/tair12/diagnosis/ADR1L1_100kb/ADR1L1_100kb.fa`
  - `benchmark/outputs/tair12/diagnosis/ADR1L1_100kb/FastNLR_config_summary.csv`
  - `benchmark/outputs/tair12/diagnosis/ADR1L1_100kb/FastNLR_target_motifs.csv`

### TAIR12 RNL-aware default 补丁

- **规则修改**：默认 motif profile 从 classic 20 扩展为 22 个 profile，仅新增经过 ADR1、
  ADR1-L2、NRG1.1、NRG1.2 和 DAR5 验证的 motif 21、23；未加入主要产生短片段的
  motif 22、24，也未加入 TIR 扩展 25–28。
- **播种修改**：新增 `[21,1]`、`[23,1]` 两个 CC→P-loop seed 及对应 signature。
  经典 NB-ARC seed 仍保持整串精确匹配；CC seed 使用连续子串匹配，因而可命中
  `[21,23,1,...]` 内部。全局 `--motif-accept-p` 保持 1e-5。
- **精度保护**：合并/延伸后仍只有两个 motif、且无 LRR 的 `CC+P-loop` 残余片段不输出；
  含 LRR 或第三结构 motif 的候选保留。
- **单元与集成测试**：`cargo test --workspace --locked` 全部通过；新增 RNL seed、
  CC 分类、完整 CC-NBARC-LRR、无 P-loop 不播种、两 motif 片段抑制和 DAR5 三 motif
  部分位点回归测试。
- **AT4G33300 局部验证**：100 kb 局部输入由 0 个位点恢复为 1 个完整
  `CC-NBARC-LRR`，motif 顺序 `21,23,1,11x5`。
- **TAIR12 全基因组回归**：位点 213→215；未匹配位点 15→15；金标准基因
  194/198→196/198；位点精确率 92.96%→93.02%；拆分基因 12→12；融合位点 10→10；
  complete NLR 145→146。
- **RNL 家族变化**：AT4G33300 新增正确位点；ADR1、ADR1-L2、NRG1.1、NRG1.2
  的基因覆盖率提高；ADR1-L3 保持 `NBARC-LRR` 不变。
- **产物**：
  - `benchmark/scripts/S17-Evaluate-TAIR12-RNL-Patch-V1.R`
  - `benchmark/reports/TAIR12-RNL-default-patch-V1.md`
  - `benchmark/metrics/tair12_rnl_patch_summary.csv`
  - `benchmark/metrics/tair12_rnl_patch_locus_changes.csv`
  - `benchmark/metrics/tair12_rnl_patch_rnl_genes.csv`
  - `benchmark/outputs/tair12/fastnlr_dna_rnl_patch/`

### FastNLR RNL-aware 新版本综合评估

- **版本标签**：`v1.5.4-rnlpatch-20261021`。所有新评估输出使用独立版本目录和
  `new_version_*` 指标名，未覆盖既有 S3/S6/S9/S17 baseline 文件。
- **RefPlantNLR 400 条阳性 CDS**：补丁前后均检出 399/400；RNL 召回 8/8；
  架构一致率 96.49%，N 端一致率 95.99%，LRR 灵敏度 99.75%，过度切分 0，
  incomplete 7，假阴性 1。补丁前后无检测状态或架构变化。
- **评估中发现的回归及修复**：DSC1 曾因 TIR 块后的单条 motif_21 交叉命中被误判为
  `CC-NBARC-LRR`。N 端规则已改为：CC/TIR 同时出现时优先更大的连续结构块，规模相同
  再取最靠近 NB-ARC 者。DSC1 恢复为 `TIR-NBARC-LRR`，RefPlantNLR 指标回到补丁前水平。
- **TAIR12**：215 个位点，200 个匹配、15 个未匹配；检出 196/198 个金标准基因；
  recall 98.99%，locus precision 93.02%，F1 95.91%；假阳性 15、假阴性 2；
  split genes 12、merged loci 10、complete NLR 146。
- **假阴性**：AT4G23440（TN17，NB-only）与 AT5G45090（PP2-A7，TIR-only）。
- **补丁差异位点**：新增 ADR1-L1、DAR5；修改 ADR1、ADR1-L2、NRG1.1、NRG1.2
  四个既有 helper NLR 位点。6 个差异位点均与金标准基因重叠；新增位点均具 RNL CC、
  P-loop 和相应 HMM 证据。
- **版本化产物**：
  - `benchmark/scripts/S18-Evaluate-TAIR12-RNL-Version-V1.R`
  - `benchmark/reports/FastNLR-RNL-version-evaluation-V1.md`
  - `benchmark/metrics/new_version_refplantnlr_*_v1.5.4-rnlpatch-20261021.csv`
  - `benchmark/metrics/new_version_tair12_*_v1.5.4-rnlpatch-20261021.csv`
  - `benchmark/outputs/refplantnlr_accuracy/patched_default_v1_5_4_rnl/`
  - `benchmark/outputs/tair12/fastnlr_dna_rnl_patch/`

### ADR1/RNL 家族分类修正（2026-10-22）

- **分类口径**：将 TAIR12 金标准从原先的“NB-ARC + CC 一律 CNL”修正为：
  TNL = NB-ARC + TIR；RNL = NB-ARC + RPW8(PF05659)；CNL = NB-ARC + 非 RPW8 CC。
- **金标准重算**：198 个金标准基因保持总数不变，分类更新为 TNL 98、RNL 6、CNL 21、
  NB-only 41、TIR-only 32。ADR1、ADR1-L1、ADR1-L2、DAR5、NRG1.1、NRG1.2 均改为 RNL。
- **软件输出**：蛋白模式 `NlrArchitecture` 新增 `family` 字段，并追加到 `.nlr.tsv` 末列。
  ADR1 保持结构类别 `CC-NBARC-LRR`，但 `family=RNL`；RPM1 与 Prf 为 `family=CNL`。
- **评估逻辑**：S6 通过位点内的 motif_21/23 证据区分 DNA 模式的 RNL 与 CNL；S7/S8
  的生物学分组增加 RNL。S5/S6/S7/S8 已重新生成 canonical 指标与图件。
- **新版本复评**：S18 使用 `v1.5.4-rnlpatch-20261022` 独立输出；TAIR12 recall 98.99%、
  locus precision 93.02%、F1 95.91%、RNL 召回 6/6、类型一致率 87.76%；
  FP 15、FN 2，split genes 12、merged loci 10。
- **报告**：`benchmark/reports/FastNLR-RNL-version-evaluation-V2.md`。

### FastNLR 1.6.0 运行时错误修复与回归（2026-10-22）

- **checkpoint 一致性**：新增版本化 `manifest.json`（format v2），记录 FastNLR 版本、输入与
  mot/store SHA-256、阈值、切片、全部组装参数、柔性/宽松播种及外部 motif 规则。恢复时使用
  当前已验证的 `AssembleParams`；缺失、损坏、未完成、输入或配置不一致均直接报错。
- **确定性修复**：checkpoint 固化首次运行时间戳，恢复时沿用；解决同一配置首次与恢复运行
  在位点起点/分类及 GFF `##date` 上的差异。TAIR12 `--no-flexible-seed --relaxed-seed 3`
  实跑首次与恢复均为 172 个位点，`.nlr.txt`、`.nlr.gff`、`.nlr.bed` 逐字节一致。
- **算法修复**：单 motif 可在追加第二个 motif 前成为 seed，`--flexible-seed 1` 生效；
  修复反向链串联尾段吸收方向；补充 Motif 与 MotifList 的 frame/ID/end/position 等确定性
  排序键；短序列按每个 motif 的实际长度独立扫描。
- **RNL/CNL 判定**：蛋白模式仅当真正被选中的 N 端 CC 含 PF05659 时判为 RNL，下游 RPW8
  交叉命中不再把 Rx_N/R1 型 CNL 误判。ADR1 为 `class=CC-NBARC-LRR, family=RNL`，RPM1 为
  `class=CC-NBARC-LRR, family=CNL`；summary 与 stats 增加 family 汇总。
- **严格错误处理**：FASTA/chopper 传播 gzip、截断和 I/O 错误；拒绝非法窗口、线程、batch、
  p 值、阈值倒置、空/非法/重复 seed/signature/category、未知 motif ID、非法 mot/store
  数字、重复键及非单调 CDF；图件和依赖基因组的输出失败不再被吞掉。
- **测试**：`cargo test --workspace --locked`、Clippy `-D warnings`、`cargo fmt --check` 和
  `git diff --check` 全部通过；release 构建成功，`fastnlr --version` 为 `1.6.0`。
- **RefPlantNLR 回归**：输出至
  `benchmark/outputs/refplantnlr_accuracy/patched_default_v1_6_0_rnl/`；检出 399/400，
  RNL 8/8，架构一致率 96.49%，过度切分 0，incomplete 7，与 1.5.4 RNL 补丁基线一致。
- **TAIR12 回归**：输出至 `benchmark/outputs/tair12/fastnlr_dna_rnl_patch_v1_6_0/`；
  215 个位点、匹配 200、未匹配 15、检出 196/198、locus precision 93.02%、
  拆分基因 12、融合位点 10、complete 146；6 个金标准 RNL 全部检出。
- **范围**：未修改 `benchmark/scripts/`，未覆盖既有评估输出。conda 版本与 SHA 待正式
  `v1.6.0` tag 生成后更新。
