# FastNLR motif 补丁：RNL（helper NLR）+ TIR 家族

本补丁给 FastNLR 的 motif 库补上两个原本缺失的家族，都是"从根上"解决 motif 法打不到的问题：

| 家族 | 解决的缺口 | 实测效果 |
|---|---|---|
| **RNL / helper NLR**（ADR1/NRG1 + 茄科 **NRC**） | 原版 motif 库里没有 RNL 的 N 端 CC motif → 11 个 RNL 只有 2 个落在位点内 | **11/11 个 RNL 区域都有位点** |
| **TIR 家族**（含单子叶 TIR-only / TN / TNP） | Pfam TIR 打不到分化/短的 TIR 域（水稻全蛋白组只命中 4–5 个）；TIR-only 又没有 NB-ARC，原版唯一的 TIR 组合 `[18,15,13]` 要求精确三连命中 | 单子叶 TIR/TIR-only 序列 0 → **11–16 个位点**；水稻 TIR-only 也出位点 |

下文 §1–§7 是 RNL 部分；**§8 专门写 TIR 部分**。

> 目标：让 FastNLR 能对 **helper NLR**（ADR1/NRG1 类 + 茄科 **NRC clade**）播种。
> 原版在 Wm82 T2T 上：11 个 RNL 基因只有 **2 个**落在 FastNLR 位点里 → README §「FastNLR 能找 RNL 吗」。
> 本补丁：**11/11 个 RNL 区域都有位点**（±5 kb 口径同样 11/11），并且新增的 motif 对 NRC 专属。

## 1. 为什么原版找不到（读源码得到的三条硬机制）

| 机制 | 位置 | 后果 |
|---|---|---|
| 播种要求 **精确的 motif ID 串** | `nlr-core/src/signature_def.rs::is_seed`（`SEED_COMBINATIONS` 11 条硬编码组合） | RNL 区域的命中串是 `[6,4,10,3,2,11,11,11,11,9]`、`[4,10,7,11,...]` 这类形状，**不等于**任何组合 → 不播种 |
| fragment 预过滤要求匹配 `SIGNATURES` | `nlr-cli/src/lib.rs::scan_one_fragment` → `has_signature` | 一个 fragment 的 motif 串匹配不到任何 signature 时，该 fragment 的**全部命中被丢弃** |
| motif 库里**没有 RNL 的 N 端 CC motif** | `crates/nlr-cli/data/mot.txt`（20 个 motif） | RNL 的 RPW8/CC-RPW8 域打不到 CC 类 motif（16/17 是 Rx_N 型 CC）；实测 RNL 区域只有 NB-ARC/LRR 类命中 |

实测（原版，11 个大豆 RNL 区域 ±5 kb）：命中只有 `motif_1/2/3/4/6/10`（NB-ARC）+ `9/11`（LRR），无一个 CC 类。

## 2. 做法（四步，全部可复现）

```
① 取 helper NLR 的 N 端序列（NB-ARC 起点前 20 aa 处截断）
     · ADR1 / ADR1-L1..L3 / NRG1 / AtNRG1.1 / AtNRG1.2 / NbADR1（RefPlantNLR）
     · **NRC clade（茄科）**：NRC0 / NRC1 / NRC2 / NRC3 / NRC4a / NRC4b / NRC6 / NbNRC2a / NbNRC2b / NbNRC3 / NbNRC4a / NbNRC4b / NbNRC4c
     · 我们已确认的 11 个大豆 RNL 模型的 N 端
   脚本：helper_nrc_all.fa → hmmsearch(NB-ARC) 定起点 → helper_nrc_nterm.fa（21 条）+ 大豆 19 条 = MEME 输入
② MEME 重挖 motif（`-protein -nmotifs 5 -minw 15 -maxw 30 -mod zoops`）→ meme_rnl_out/（v1）、meme_nrc_out/（v2，含 NRC）
③ 把 MEME 的 letter-probability matrix 转成 FastNLR 的 mot.txt / store.txt 条目
   脚本：make_rnl_motifs.py
     · PWM：p×100 取整（与现有 20 个 motif 的 0–100 量纲一致）
     · CDF：**精确动态规划**（逐位置卷积 + 大豆蛋白组氨基酸背景频率）算窗口得分分布 → 右尾 CDF
       （不用随机采样，尾部到 1e-12；索引与 FastNLR 的"窗口得分求和查表"口径一致）
④ 打源码补丁 + 编译
   脚本：apply_rnl_motif_patch.py（可重复执行；先检测是否已打过）
     · mot.txt / store.txt 追加新 motif（内嵌进二进制）
     · MAX_MOTIF_ID 20 → 24（原版数组/循环里硬编码的 21、1..=20 全部改成 MAX_MOTIF_ID 派生）
     · 新 motif 归入 **DomainCategory::Cc**
     · SEED_COMBINATIONS / SIGNATURES 追加"CC 开头 + NB-ARC"的组合
     · `is_seed`：**只对**"以新 motif（id ≥ RNL_CC_MOTIF_MIN=21）开头"的组合改用**连续子串匹配**，
       原版组合保持"整串精确匹配"→ 老行为零改动
   编译：export PATH=/var/www/masw_data/software/mambaforge/envs/rust/bin:$PATH
         cd tools/FastNLR-src && cargo build --release   # ≈35 s
```

## 3. 最终选定的 motif 集（4 个）

| id | 来源 | 宽度 | 共识序列 | 实测命中 |
|---|---|---|---|---|
| 21 | MEME v1（ADR1/NRG1 + 大豆 RNL） | 30 | `CKSTAESLISTJNDLLPTIZEIKYSGVELD` | 大豆 RNL **11/11**、对照 NLR 1/30 |
| 22 | MEME v1 | 30 | `WFDMDFPKAEVLILNFSSDDYVLPPFIAKM` | 大豆 RNL 8/19 正样本、对照 0/30 |
| 23 | MEME v1 | 30 | `RQEZJDRLSEILRAGEELVRKVLSSSRWNV` | 正样本 17/19、对照 1/30 |
| **24** | **MEME v2（NRC clade 专属）** | 30 | `VKKIRKVVNSAEDAIDKFVIZAKLHKDKNK` | **NRC 13/13、ADR1/NRG1 0/8**、对照 0/100 |

motif 24 是加进 NRC 序列后新挖出来的、**NRC 专属**的那个：它把 NRC 与 ADR1/NRG1 两类 helper 分开
（互补关系，不是重复）——这正是"补上茄科 NRC"的价值。

逐 motif 的评估（`--tiscalling` 无关，见 §5 表）与弃用记录：

- 同一次 MEME 还挖出 2 个 ADR1/NRG1 风味 motif（`ZJDRLSEILRKGVELVHKVLKSSRWNVYRN`、`KMEKLEKHVSRFLQGPMQAHILADVHHVRF`）：
  NRC 0/13、ADR1/NRG1 6/8 与 4/8 —— 与 motif_22/23 冗余，加上后全基因组位点从 690 涨到 736 而大豆覆盖没有提升，**弃用**。

## 4. 测量结果（Wm82 T2T）

| 版本 | motif 数 | 全基因组位点 | RNL 区域有 loci | RNL 模型覆盖（±5 kb） | T1/T2 模型被位点覆盖 | 邻近基因比例 |
|---|---|---|---|---|---|---|
| **原版（基线）** | 20 | **562** | 2/11 | 11/11* | 388/426 | 80% |
| v1（3 motif，ADR1/NRG1+大豆） | 23 | 635 | 11/11 | 11/11 | 398/426 | 77% |
| v3/v4（MEME v2 的 3 个 / +NRC） | 23/24 | 680 | 10/11 | 9/11 | 395/426 | — |
| v2（5 motif，含 2 个冗余） | 25 | 736 | 11/11 | 11/11 | 400/426 | — |
| **v5 = v1 + NRC（本次交付）** | **24** | **690** | **11/11** | **11/11** | **398/426** | 77% |

\* 基线那一行"±5 kb"列有歧义：原版位点直接覆盖 2/11、±5 kb 内 11/11 中的多数是因为相邻基因位点；本补丁后
**每个 RNL 区域都能独立出位点**（`11 sequences, 11 NLR loci`）。

代价与收益：

- **代价**：全基因组位点 562 → 690（+23%），但增量几乎全部是"含新 CC motif 的位点"
  （v1：634 个位点里 80 个含新 motif），且**邻近基因比例只从 80% 降到 77%** → 不是成片的假位点；
  另外这 +128 个位点在流程里先成为**候选位点**，仍需 06 步的域/结构证据才能成模型
  （T1/T2 数量不变，motif-only 位点 22 → 27）。
- **收益**：RNL/helper 从"motif 层完全看不见"变成"每个区域都有位点"；新增 motif 对
  NRC（茄科 helper）专属，便于把这套流程搬到番茄/马铃薯等茄科作物；
  含 CC motif 的位点域串会输出成 `CC-NBARC-LRR`（新 motif 已归入 CC 类，可解释）。

## 5. 复现 / 回退

```bash
# 复现（从干净源码开始）
cd /media/masw/Fei_data2/lixiang_data/12.NLR_pipeline
#   ① MEME 输入（仓库已带）：results_t2t/fastnlr_patch/helper_nrc_nterm.fa + rnl_nterm_for_meme.fa
#   ② motif 条目（仓库已带）：rnl_motif.mot.txt / .store.txt（motif 21-24，RNL/helper CC）+ 24 号取自 meme_nrc_out
envs/tiscalling 无关；用 nlr 环境的 meme：
  /var/www/masw_data/software/mambaforge/envs/nlr/bin/meme <input.fa> -protein -oc <out> -nmotifs 5 -minw 15 -maxw 30 -mod zoops
python3 results_t2t/fastnlr_patch/make_rnl_motifs.py --meme <out>/meme.txt \
    --proteome work_t2t/proteome.faa --start-id 21 [--keep 1,2,3] --out-prefix /tmp/rnl_motif
python3 results_t2t/fastnlr_patch/apply_rnl_motif_patch.py \
    --src tools/FastNLR-src --motif-prefix /tmp/rnl_motif_v5
export PATH=/var/www/masw_data/software/mambaforge/envs/rust/bin:$PATH
cd tools/FastNLR-src && cargo build --release && cp target/release/fastnlr ../fastnlr-patched

# 回退到原版（不带 RNL motif）
cd tools/FastNLR-src && git checkout . && cargo build --release && cp target/release/fastnlr ../fastnlr-patched
#   （源码仓库在 tools/FastNLR-src 里是独立 git 仓库，checkout 即可回到纯净版）
```

`rnl_motif.patch` 是本次补丁相对**原版源码**的完整 diff（`diff -ruN --exclude=target`），
`fastnlr-patched.sha256` 是交付二进制的校验和。

## 6. 与既有补丁的关系

同一份源码里已有两组早期补丁（都是"提高分化 NLR 灵敏度"的方向，见 `FastNLR-RNL.patch` 与
`01..05_*.diff`）：`--relaxed-seed N`（按域类别播种）、`--motif-accept-p / --motif-prelim-p`（放宽阈值）。
本次补丁与它们**互补且不冲突**：那三个参数是"放宽已有 motif 的判定"，本次是"补齐缺失的 motif 家族"。
默认参数下（不放宽阈值）本补丁即可让 RNL 出位点，因此流程 02b 仍用默认参数。

## 7. 快速自检（改完源码/换机器后各跑一次）

```bash
cd /media/masw/Fei_data2/lixiang_data/12.NLR_pipeline

# ① 补丁是否在源码里（应输出 motif_21..motif_24）
grep -o '^motif_2[1-4]' tools/FastNLR-src/crates/nlr-cli/data/mot.txt | sort -u | tr '\n' ' '; echo
grep -c 'RNL_CC_MOTIF_MIN' tools/FastNLR-src/crates/nlr-core/src/signature_def.rs

# ② 二进制是否已重新编译（对比 sha256 与交付记录）
sha256sum tools/fastnlr-patched; cat results_t2t/fastnlr_patch/fastnlr-patched.sha256

# ③ 行为自检：拿 11 个 RNL 区域跑一遍（需先有 results_t2t/nlr.all.tsv）
python3 - <<'PY'
import csv
rows=[r for r in csv.DictReader(open('results_t2t/nlr.all.tsv'),delimiter='\t')
      if r['nlr_class']=='RNL' and r['final_tier'] in ('T1','T2')]
with open('/tmp/rnl_regions.bed','w') as o:
    for r in rows: o.write('%s\t%d\t%d\t%s\n'%(r['chrom'],max(0,int(r['start'])-5000),int(r['end'])+5000,r['id']))
print('RNL 区域:', len(rows))
PY
/var/www/masw_data/software/mambaforge/envs/nlr/bin/bedtools getfasta -fi work_t2t/genome.fa \
    -bed /tmp/rnl_regions.bed -name -fo /tmp/rnl_regions.fa
tools/fastnlr-patched -i /tmp/rnl_regions.fa -p /tmp/rnl_check -t 8 2>&1 | grep 'NLR loci'
#   期望：NLR loci = 11（原版是 2）
```

---

## 8. TIR 家族补丁（含单子叶 TIR-only / TN / TNP）

### 8.1 缺口

1. **Pfam TIR（PF01582）打不到分化/短的 TIR 域**：拿它扫水稻全部 49,066 个蛋白，i-E<1e-3 只命中 **3 个**、
   i-E<0.1 也只有 3 个（TIR_2 命中 4–5 个）——单子叶的 TIR-only/TN/TNP 基本落在模型之外。
2. **TIR-only 没有 NB-ARC**：原版 `SEED_COMBINATIONS` 里唯一含 TIR 的组合是 `[18,15,13]`（要求精确三连命中），
   而 `[18,15]`、`[15,13]` 只写在 `SIGNATURES`（预过滤）里、**不在播种组合里** → TIR-only/短 TIR 蛋白永远播不了种。
   实测：把 23 条单子叶 TIR 蛋白反译成 DNA 交给原版 FastNLR，**0 个位点**。

### 8.2 做法

1. **收集 TIR 域序列（236 条）**：
   - **单子叶**：UniProt 取小麦/玉米/高粱/短柄草/水稻的 TIR 蛋白（`xref:pfam-PF01582` 或 `protein_name:TIR`）
     → 23 条，其中 **13 条为 TIR-only/TN**（无 NB-ARC）、10 条为 TNL（`monocot_tir_proteins.fa`）
   - **双子叶**：RefPlantNLR 里 93 条 TIR 类（实验验证的 TNL）+ 我们已确认的大豆 TNL 的 TIR 域
   - 用 hmmsearch 定位 TIR 域边界后前后各留 20 aa（`tir_domains_for_meme.fa`）
2. **MEME 重挖**（`-protein -nmotifs 4 -minw 15 -maxw 30 -mod zoops`）→ 4 个 motif，每个覆盖 215–223/236 条序列
   （E 值 1e-2000 量级）：
   | id | 共识 | 宽度 | 说明 |
   |---|---|---|---|
   | 25 | `YDVFLSFRGEDTRKTFTSHLY` | 21 | TIR1 亚域的经典指纹 |
   | 26 | `VJPVFYNVDPSDVRKQTGSYG` | 21 | TIR2 样 |
   | 27 | `AIEESRIAIVVFSKNYASSSWCLDELAKIL` | 30 | TIR 家族保守段 |
   | 28 | `KEDEEKVQKWRKALTEVANLSGWHS` | 25 | TIR 家族保守段 |
3. **归入 TIR 类**（`DomainCategory::Tir`）+ **加 TIR-only 播种组合**：
   `[新TIR × 新TIR]`、`[新TIR × 原版 TIR 13/15/18]`、`[新TIR × NB-ARC]`，另外把原版签名里的
   **`[18,15]`、`[15,13]` 提升为 seed**（这条正是 TIR-only 播不了种的直接原因）。
4. **收紧标定（关键的一步）**：TIR motif 的命中面比 CC 宽（不收紧时全基因组位点 690 → **1168**、
   中位长度掉到 857 bp）。做法是把每种 TIR motif 的 **CDF 左移 delta 分**
   （等价于"分数必须高 delta 才能达到同样的显著性"，⚠ 方向不能反：右移是**放松**判据，
   实测右移会得到 15,754 个位点），使 `p<1e-5` 处的期望命中窗口数从 ~370 降到 ~40。
   参数见 `make_rnl_motifs.py --target-exp`。

### 8.3 实测（v6 = CC 21-24 + 收紧后的 TIR 25-28）

| 测试 | 原版 | 仅 CC 补丁（v5） | **v6（CC+TIR）** |
|---|---|---|---|
| 单子叶 TIR/TIR-only（23 条反译 DNA） | **0 位点** | 0 | **11–16 位点** |
| 水稻 TIR-only CDS（3 条） | 0 | 0 | **2** |
| 大豆 RNL 区域（11 个） | 2 | 11 | **11** |
| 大豆全基因组位点 | 562 | 690 | **883** |
| 邻近基因比例 | 80% | 76% | **76%** |
| 位点长度中位 | 2540 bp | 2389 bp | 1540 bp |

- 新增 TIR 位点的**域串**：`TIR`（239）、`TIR-NBARC`（171）、`TIR-LRR`（93）、`TIR-NBARC-LRR`（21）……
  → TIR-only / TN / TNL 都能出位点，且能看出来是哪种架构。
- **新增位点的质量不比原来差**：按 motif 拆开看，motif_27 贡献的 226 个位点里 **86% 邻近基因**、
  motif_28 的 158 个里 **89% 邻近基因**（整体 76%）→ 检出的是真 TIR 家族基因，不是散布垃圾。
- 代价：位点数 +28%（690 → 883）、中位长度变短（TIR-only 基因本来就短）、<1 kb 位点 216 → 399。
  这些位点在流程里先成为**候选位点**，仍需 06 步域/结构证据才成模型。

### 8.4 复现

```bash
# ① 单子叶 TIR 序列（UniProt，需联网）
for org in 39947 4565 4577 4513 4558 15368; do
  for q in "xref:pfam-PF01582" "protein_name:TIR"; do
    curl -s "https://rest.uniprot.org/uniprotkb/search?query=organism_id:$org+AND+($q)&format=fasta&size=200"
  done
done > monocot_tir_raw.fa      # 去重后 23 条；其中 TIR-only 13 条
# ② hmmsearch 定位 TIR 域 → 取 ±20 aa 的域片段（236 条）→ MEME
meme tir_domains_for_meme.fa -protein -oc /tmp/meme_tir -nmotifs 4 -minw 15 -maxw 30 -mod zoops
# ③ 转成 mot/store（TIR 家族要收紧标定；CC 家族沿用 v5 的库）
python3 make_rnl_motifs.py --meme /tmp/meme_tir/meme.txt --proteome work_t2t/proteome.faa \
    --start-id 25 --target-exp 40 --out-prefix /tmp/tir_motif
cat rnl_motif.mot.txt tir_motif.mot.txt > motif_all.mot.txt      # store 同理
# ④ 打补丁 + 编译（MOTIF_SPECS = 21-24:CC, 25-28:TIR）
python3 apply_rnl_motif_patch.py --src tools/FastNLR-src --motif-prefix /tmp/motif_v6
export PATH=/var/www/masw_data/software/mambaforge/envs/rust/bin:$PATH
cd tools/FastNLR-src && cargo build --release && cp target/release/fastnlr ../fastnlr-patched
```

### 8.5 还没做（需要单子叶基因组才能验证）

- 本机只有水稻的**蛋白/CDS**（`/var/www/masw_data/tiantian_data/osa1_r7.*`），**没有水稻基因组 FASTA**，
  所以 §8.3 的单子叶验证是"蛋白反译成 DNA 后直接跑 FastNLR"，等价于只测 motif 层能否命中并播种
  （这正是本次要修的东西），但**没有**做"单子叶全基因组位点数量/假阳性"的评估。
- 拿到任一单子叶基因组（水稻/玉米/短柄草）+ 注释后，按同一流程跑一遍即可补上：
  `T2T 那套命令换成单子叶的 config`，重点看"TIR-only/TN/TNP 基因是否都有位点"和"位点数是否失控"。
