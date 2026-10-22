# FastNLR 柔性播种升级方案（category / rank-based seeding）

> 作者：Jiwen Zhao (https://github.com/CropCoder)
> 日期：2026-09-23
> 状态：已实现；v1.5.x 中默认开启
> 关联：RefPlantNLR 阳性集评测（`benchmark/metrics/refplantnlr_*.csv`）

> **现行行为（2026-10-15 后）**：默认启用 category/rank 柔性播种，
> `min_nbarc = 3`、`require_ploop = false`，并要求 NB-ARC rank **严格递增**。
> 柔性条件与精确组合取并集；`--no-flexible-seed` 可回退到仅精确组合。
> 下文保留设计演进和修订记录，旧的“单调不减”代码仅用于解释历史问题。

## 1. 背景与问题

柔性播种加入前，`AnnotatorSignatureDefinition::is_seed` 对内置 motif（1..=20）采用
**整串精确相等**（`s == seq`）。`find_seeds` 从某个 motif 起，只要相邻 motif 的
`dna_start` 差值 ≤500 bp 就持续 append 进 `s`，每 append 一次查一次 `is_seed(s)`，命中即播种。

结果是：真实 NLR 的 motif 命中串只要与 11 条内置组合不完全一致——少一个、多一个、
换顺序——就永远无法播种。典型反例是 RNL（helper NLR），其命中串形如
`[6,4,10,3,2,11,11,9]`、`[4,10,7,11,...]`，不等于任何内置组合。

当时的实测后果（RefPlantNLR 400 条阳性 CDS，历史数据）：

| 配置 | 基因级召回 | RNL 召回（共 8 个） |
|------|-----------|---------------------|
| default（精确匹配） | 393/400（98.25%） | 3/8 |
| extended（外部库加组合补丁） | 398/400（99.50%） | 8/8 |

扩展库靠「继续往 SEED_COMBINATIONS / SIGNATURES 里加组合」打补丁，本质是在补
穷举的窟窿，没有改掉「精确匹配」这个前提。本方案从根上解决。

## 2. 核心思路

把匹配对象从「具体 motif ID」换成「结构域类别 + rank 顺序」。

内置 11 条种子组合本质上是 NB-ARC 域内 motif 按 rank 序（`[1,6,4,5,10,3,12,2]`）的
连续片段。与其穷举这些片段，不如直接表达这条规律：

> 一个片段里，若出现的 NB-ARC 类 motif 数量 ≥ min_nbarc，且它们的 rank 严格递增
> （按 P-loop → 末端的正确顺序），即视为种子。

RNL 的 `[6,4,10,3,2]` 的 rank 是 `5,6,8,9,11`，严格递增，因此能播种。

## 3. 算法设计

### 3.1 配置结构

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FlexibleSeedConfig {
    pub min_nbarc: usize,      // 最少 NB-ARC motif 数，默认 3
    pub require_ploop: bool,   // 是否要求含 P-loop（motif_1），默认 false
}

impl Default for FlexibleSeedConfig {
    fn default() -> Self {
        FlexibleSeedConfig { min_nbarc: 3, require_ploop: false }
    }
}
```

> `require_ploop` 默认 false 的原因：目标 RNL 命中串常不含 motif_1，若强制 P-loop
> 会导致「救了等于没救」。它作为可选的精度旋钮保留。

### 3.2 柔性播种判定

```rust
pub fn is_seed_flexible(&self, seq: &[MotifId], cfg: &FlexibleSeedConfig) -> bool {
    // 只取 NB-ARC 类 motif，保持出现顺序
    let nbarc: Vec<MotifId> = seq.iter().copied()
        .filter(|&id| self.is_nbarc(id))
        .collect();
    if nbarc.len() < cfg.min_nbarc { return false; }
    if cfg.require_ploop && !nbarc.contains(&PLOOP_MOTIF) { return false; }
    // rank 严格递增
    nbarc.windows(2).all(|w| self.rank(w[0]) < self.rank(w[1]))
}
```

### 3.3 柔性预过滤

```rust
pub fn has_signature_flexible(&self, ids: &[MotifId], cfg: &FlexibleSeedConfig) -> bool {
    // 预过滤面对整个片段，保留宽松判据，避免尾部重复 motif 丢弃合法种子前缀
    self.has_signature(ids) || self.is_seed_flexible_loose(ids, cfg)
}
```

预过滤必须一起放宽，否则 RNL 片段在扫描阶段就被 `has_signature` 丢弃，根本到不了
播种这一步（这是 RNL 漏检的第一道闸）。

## 4. 改动清单

| 文件 | 改动 |
|------|------|
| `crates/nlr-core/src/signature_def.rs` | 新增 `FlexibleSeedConfig`；`AnnotatorSignatureDefinition` 加字段 `flexible_seed: Option<FlexibleSeedConfig>`、`with_flexible_seed()`、`flexible_seed()`、`is_seed_flexible()`、`has_signature_flexible()` |
| `crates/nlr-assemble/src/lib.rs` | `find_seeds` 播种判定改为三态：`flexible_seed()` 时取 `is_seed_flexible() || is_seed()`（并集）→ 否则 `relaxed_seed_min_nbarc` → 否则 `is_seed` |
| `crates/nlr-cli/src/lib.rs` | `RunConfig` 加 `flexible_seed: Option<FlexibleSeedConfig>`；`filter_motifs_by_signature` 按需改用 `has_signature_flexible`；构造 `signature_def` 时注入 |
| `crates/nlr-cli/src/main.rs` | 新增 `--flexible-seed <n>` 与 `--require-ploop`，写入 `config.flexible_seed` |

`find_seeds` 播种判定伪代码：

```rust
let seeded = if let Some(cfg) = def.flexible_seed() {
    // 并集：柔性条件与精确组合同时生效（见下方修订说明）
    def.is_seed_flexible(&s, &cfg) || def.is_seed(&s)
} else if let Some(n) = relaxed_seed_min_nbarc {
    def.is_seed_relaxed(&s, n)
} else {
    def.is_seed(&s)
};
```

> **修订（2026-10-12）**：最初实现为 `flexible_seed()` 时**只**看柔性条件（替换语义），
> 与 3.3 节的预过滤（`has_signature || is_seed_flexible`，并集）不一致：命中精确签名而被
> 保留的 motif 簇，在播种阶段却不认精确组合，于是「保留下来却组装不出位点」。实测在
> 142 Mb 的真实基因组上，替换语义会丢掉默认配置能检出的 38/207 个位点（-18%），
> 且 RefPlantNLR 的 Pb1（精确组合 `[12,2,8]`，NB-ARC motif 仅 2 个）在
> `--flexible-seed 3` 下反而漏检。改为并集后：38 个位点全部找回，金标准召回
> 99.50% → 99.75%，过度切分不变（3 个基因）。并集只会新增种子，默认路径行为逐字节不变。

> **修订（2026-10-15）**：播种判据由「rank **单调不减**」收紧为「rank **严格递增**」
> （`<=` → `<`）。8 个 NB-ARC motif 的 rank 两两不同，因此两种判据只在"同一个 motif 被
> 重复命中"时才有差别——典型是低复杂度区域里连续三次 `motif_2`（rank 11,11,11）。
> 这类重复不构成结构域的先后顺序，实测会在基因间区造出假位点（TAIR12 上的
> `OZ408684.1_nlr2`，619 bp、无 P-loop、无 LRR），且在与 NLR-Annotator 的横向对比里
> 是唯一可归因于实现层面的假阳性。收紧后 TAIR12 上位点 214→213、召回不变（194/198）、
> 位点级假阳性率 7.48%→7.04%，RefPlantNLR 四档指标逐字节不变。
>
> **注意**：片段预过滤 `has_signature_flexible` 仍用宽松判据（`is_seed_flexible_loose`）。
> 预过滤面对的是**整个片段的命中串**，若一并收紧，片段尾部的一个重复 motif 会让该片段
> 里合法的种子前缀被整段丢弃——实测会误删 TAIR12 上 `motif 1,4,10,1` 的 CC-NBARC-LRR 位点。

## 5. CLI 接口

```bash
# 默认配置：至少 3 个 NB-ARC motif 按 rank 严格递增即播种
fastnlr -i genome.fa -o out

# 提高精度：额外要求含 P-loop
fastnlr -i genome.fa -o out --require-ploop

# 提高召回：降到 2 个 NB-ARC motif
fastnlr -i genome.fa -o out --flexible-seed 2
```

## 6. 可靠性设计

三档谱系与精度/召回权衡：

| 档位 | 规则 | 召回 | 精度 |
|------|------|------|------|
| 精确 ID 匹配（关闭柔性播种时） | `s == combo` | 低 | 高 |
| 类别 + 严格递增 rank（当前默认） | NB-ARC≥min 且 rank 严格递增，可选 P-loop | 高 | 可控 |
| 只数类别（现有 `--relaxed-seed`） | NB-ARC≥n，不看顺序 | 最高 | 最低 |

可靠性保障：

1. **rank 严格递增是强约束**：随机打中的不同 NB-ARC motif 恰好按正确顺序排列概率低；
   `merge_seeds` 仍使用非递减顺序合并，并允许 LRR 重排。
2. **P-loop 可选锚点**：需要更高精度时用 `--require-ploop`。
3. **min_nbarc 可调**：3 为默认，2 提召回，4+ 提精度。
4. **保留既有机制**：柔性播种与精确组合取并集，并使用 `--no-flexible-seed` 提供回退路径。

## 7. 测试计划

单元测试（`crates/nlr-core/tests/core.rs`）：

- `is_seed_flexible(&[6,4,10,3,2], cfg)` 为 true（RNL 形状）
- `is_seed_flexible(&[6,4,10,3,2], cfg{require_ploop:true})` 为 false（缺 P-loop）
- rank 反序（如 `[2,3,12,10]` 映射后非单调）为 false
- `is_seed_flexible(&[1,6,4], cfg)` 为 true（兼容原精确组合）
- `has_signature_flexible` 对 RNL 串返回 true

回归/验证（`benchmark/scripts/S3-Run-RefPlantNLR-Accuracy-V1.R`）：

- 阳性集：400 条 RefPlantNLR CDS，比较严格精确基线（`--no-flexible-seed`）、当前柔性默认和
  extended 三档的基因级召回、RNL 召回、结构域一致性、过度切分。当前 S3 中的 `dna_default`
  与 `dna_flexible` 均启用柔性播种，结果相同，解释时不能把它们当作严格基线与柔性档的重叠。
- 负对照：对 400 条 CDS 做打乱（保持碱基组成、破坏 NLR 信号），统计柔性播种引入的
  假阳性位点数，作为精度指标。
- 验收重点：柔性播种不降低精确组合检出能力，并比较 RNL 召回、整体召回、架构一致性与
  负对照假阳性。

## 8. 回滚与兼容性

- 柔性播种在 CLI 中默认开启；所有字段仍通过 builder 注入，`--no-flexible-seed` 可回退到
  仅精确组合。
- `require_ploop` 与 `min_nbarc` 均为独立开关，互不影响。

## 9. 记录

- 评测产物：`benchmark/metrics/refplantnlr_gene_table.csv`、`refplantnlr_summary.csv`、
  `refplantnlr_by_type.csv`、`refplantnlr_by_domain.csv`。
- 评测图：`benchmark/figures/FigA1~FigA4-*.pdf`。

## 10. RNL-aware default 补丁（2026-10-21）

ADR1-L1（AT4G33300）局部诊断显示，default 库只留下 `motif_1 -> motif_11xN`，无法满足
3 个 NB-ARC motif 的柔性 seed；但 helper CC 与 P-loop 相邻。补丁采用两项保守扩展：

1. default profile 新增 motif 21、23，并加入 `21,1`、`23,1` 两个精确 seed 及对应 signature；
2. 新增 CC seed 使用连续子串匹配，经典 NB-ARC seed 仍要求整串精确相等；合并/延伸后仍只有
   两个 motif 且无 LRR 的 `CC+P-loop` 残余片段不输出。
3. CC/TIR 同时出现在首个 NB-ARC 前时，N 端域优先连续块规模更大的类别，规模相同再取更靠近
   NB-ARC 者；防止 TIR 基因中的单条 RNL-CC 交叉命中覆盖完整 TIR 块。

TAIR12 全基因组回归：位点 213→215，未匹配位点 15→15，金标准基因 194→196/198，
拆分基因 12→12、融合位点 10→10。AT4G33300 从未检出恢复为完整 `CC-NBARC-LRR`。
