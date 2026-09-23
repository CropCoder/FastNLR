# FastNLR 柔性播种升级方案（category / rank-based seeding）

> 作者：Jiwen Zhao (https://github.com/CropCoder)
> 日期：2026-09-23
> 状态：设计定稿，待实现
> 关联：RefPlantNLR 阳性集评测（`benchmark/metrics/goldnlr_*.csv`）

## 1. 背景与问题

当前播种判定 `AnnotatorSignatureDefinition::is_seed` 对内置 motif（1..=20）采用
**整串精确相等**（`s == seq`）。`find_seeds` 从某个 motif 起，只要相邻 motif 距离
≤500 bp 就持续 append 进 `s`，每 append 一次查一次 `is_seed(s)`，命中即播种。

结果是：真实 NLR 的 motif 命中串只要与 11 条内置组合不完全一致——少一个、多一个、
换顺序——就永远无法播种。典型反例是 RNL（helper NLR），其命中串形如
`[6,4,10,3,2,11,11,9]`、`[4,10,7,11,...]`，不等于任何内置组合。

实测后果（RefPlantNLR 400 条阳性 CDS）：

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

> 一个片段里，若出现的 NB-ARC 类 motif 数量 ≥ min_nbarc，且它们的 rank 单调不减
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
    // rank 单调不减
    nbarc.windows(2).all(|w| self.rank(w[0]) <= self.rank(w[1]))
}
```

### 3.3 柔性预过滤

```rust
pub fn has_signature_flexible(&self, ids: &[MotifId], cfg: &FlexibleSeedConfig) -> bool {
    self.has_signature(ids) || self.is_seed_flexible(ids, cfg)
}
```

预过滤必须一起放宽，否则 RNL 片段在扫描阶段就被 `has_signature` 丢弃，根本到不了
播种这一步（这是 RNL 漏检的第一道闸）。

## 4. 改动清单

| 文件 | 改动 |
|------|------|
| `crates/nlr-core/src/signature_def.rs` | 新增 `FlexibleSeedConfig`；`AnnotatorSignatureDefinition` 加字段 `flexible_seed: Option<FlexibleSeedConfig>`、`with_flexible_seed()`、`flexible_seed()`、`is_seed_flexible()`、`has_signature_flexible()` |
| `crates/nlr-assemble/src/lib.rs` | `find_seeds` 播种判定改为三态：`flexible_seed()` 优先 → `relaxed_seed_min_nbarc` → `is_seed` |
| `crates/nlr-cli/src/lib.rs` | `RunConfig` 加 `flexible_seed: Option<FlexibleSeedConfig>`；`scan_one_fragment` 预过滤按需改用 `has_signature_flexible`；构造 `signature_def` 时注入 |
| `crates/nlr-cli/src/main.rs` | 新增 `--flexible-seed <n>` 与 `--require-ploop`，写入 `config.flexible_seed` |

`find_seeds` 播种判定伪代码：

```rust
let seeded = if let Some(cfg) = def.flexible_seed() {
    def.is_seed_flexible(&s, &cfg)
} else if let Some(n) = relaxed_seed_min_nbarc {
    def.is_seed_relaxed(&s, n)
} else {
    def.is_seed(&s)
};
```

## 5. CLI 接口

```bash
# 柔性播种：连续 3 个 NB-ARC motif 按 rank 顺序即播种
fastnlr -i genome.fa -o out.txt --flexible-seed 3

# 提高精度：额外要求含 P-loop
fastnlr -i genome.fa -o out.txt --flexible-seed 3 --require-ploop

# 提高召回：降到 2 个 NB-ARC motif
fastnlr -i genome.fa -o out.txt --flexible-seed 2
```

## 6. 可靠性设计

三档谱系与精度/召回权衡：

| 档位 | 规则 | 召回 | 精度 |
|------|------|------|------|
| 精确 ID 匹配（现状默认） | `s == combo` | 低 | 高 |
| 类别 + rank 单调（本方案） | NB-ARC≥min 且 rank 递增，可选 P-loop | 高 | 可控 |
| 只数类别（现有 `--relaxed-seed`） | NB-ARC≥n，不看顺序 | 最高 | 最低 |

可靠性保障：

1. **rank 单调是强约束**：随机打中的 motif 恰好按 NB-ARC 正确顺序排列概率低；且
   `merge_seeds` 已在用 rank 单调作为合并条件，该信号在项目内是被验证过的。
2. **P-loop 可选锚点**：需要更高精度时用 `--require-ploop`。
3. **min_nbarc 可调**：3 为默认，2 提召回，4+ 提精度。
4. **保留既有机制**：默认仍走精确匹配，柔性播种先 opt-in，验证通过再考虑设为默认。

## 7. 测试计划

单元测试（`crates/nlr-core/tests/core.rs`）：

- `is_seed_flexible(&[6,4,10,3,2], cfg)` 为 true（RNL 形状）
- `is_seed_flexible(&[6,4,10,3,2], cfg{require_ploop:true})` 为 false（缺 P-loop）
- rank 反序（如 `[2,3,12,10]` 映射后非单调）为 false
- `is_seed_flexible(&[1,6,4], cfg)` 为 true（兼容原精确组合）
- `has_signature_flexible` 对 RNL 串返回 true

回归/验证（R 脚本，复用 `benchmark/Script5-*` 模式）：

- 阳性集：400 条 RefPlantNLR CDS，对比 default / extended / flexible 三档的
  基因级召回、RNL 召回、结构域一致性、过度切分。
- 负对照：对 400 条 CDS 做打乱（保持碱基组成、破坏 NLR 信号），统计柔性播种引入的
  假阳性位点数，作为精度指标。
- 目标：flexible 档 RNL 召回接近 8/8，且负对照假阳性不显著高于 default 档。

## 8. 回滚与兼容性

- 柔性播种默认关闭，所有字段用 builder 注入，旧路径（`is_seed` / `relaxed_seed`）
  行为零改动，可随时回滚。
- `require_ploop` 与 `min_nbarc` 均为独立开关，互不影响。

## 9. 记录

- 评测产物：`benchmark/metrics/goldnlr_gene_table.csv`、`goldnlr_summary.csv`、
  `goldnlr_confusion.csv`。
- 评测图：`benchmark/figures/Fig8~Fig10-goldnlr-*.pdf`。
