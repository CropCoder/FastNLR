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
