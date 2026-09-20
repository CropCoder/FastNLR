//! NLR static rule tables and signature definitions.
//!
//! The tables define:
//! - rank, category (domain category), and RGB color for each of the 20 motifs;
//! - 11 seed combinations (findSeeds targets);
//! - 18 "NLR signatures" (contiguous motif sequence patterns used for pre-filtering);
//! - consensus sequences for 8 NB-ARC motifs;
//! - P-loop motif id (motif_1).

/// Motif id（内置 1..=20；外部 mot.txt 可带更多，见 BUILTIN_MOTIF_COUNT）。
///
/// Built-in motifs use ids 1..=20. The upper bound is a constant so external `mot.txt`/`store.txt`
/// files can carry additional motifs (for example RNL/helper-CC and TIR family profiles).
pub type MotifId = u8;


/// Domain category.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DomainCategory {
    Nbarc,
    Lrr,
    Tir,
    Cc,
    Linker,
    Na,
}

impl DomainCategory {
    /// Construct from a category string.
    pub fn from_str(s: &str) -> Self {
        match s {
            "NBARC" => DomainCategory::Nbarc,
            "LRR" => DomainCategory::Lrr,
            "TIR" => DomainCategory::Tir,
            "CC" => DomainCategory::Cc,
            "LINKER" => DomainCategory::Linker,
            _ => DomainCategory::Na,
        }
    }

    /// Convert back to a category string for output.
    pub fn as_str(self) -> &'static str {
        match self {
            DomainCategory::Nbarc => "NBARC",
            DomainCategory::Lrr => "LRR",
            DomainCategory::Tir => "TIR",
            DomainCategory::Cc => "CC",
            DomainCategory::Linker => "LINKER",
            DomainCategory::Na => "NA",
        }
    }
}

/// P-loop motif id.
pub const PLOOP_MOTIF: MotifId = 1;

/// Rank table (index = motif id; `RANKS[0]` is an unused placeholder).
pub const RANKS: [u8; 29] = [
    0, 4, 11, 9, 6, 7, 5, 13, 12, 14, 8, 14, 10, 3, 3, 2, 2, 1, 1, 14, 15,
    14,     14,     14,     14,     14,     14,     14, 14,   // 8 个 RNL(helper) 专用 CC motif（不参与 NB-ARC 排序）
];

/// Category table (index = motif id).
pub const CATEGORIES: [DomainCategory; 29] = [
    DomainCategory::Na,
    DomainCategory::Nbarc,  // 1
    DomainCategory::Nbarc,  // 2
    DomainCategory::Nbarc,  // 3
    DomainCategory::Nbarc,  // 4
    DomainCategory::Nbarc,  // 5
    DomainCategory::Nbarc,  // 6
    DomainCategory::Linker, // 7
    DomainCategory::Linker, // 8
    DomainCategory::Lrr,    // 9
    DomainCategory::Nbarc,  // 10
    DomainCategory::Lrr,    // 11
    DomainCategory::Nbarc,  // 12
    DomainCategory::Tir,    // 13
    DomainCategory::Na,     // 14
    DomainCategory::Tir,    // 15
    DomainCategory::Cc,     // 16
    DomainCategory::Cc,     // 17
    DomainCategory::Tir,    // 18
    DomainCategory::Lrr,    // 19
    DomainCategory::Na,     // 20
    DomainCategory::Cc,     // 21（新增 CC 家族 motif）
    DomainCategory::Cc,     // 22（新增 CC 家族 motif）
    DomainCategory::Cc,     // 23（新增 CC 家族 motif）
    DomainCategory::Cc,     // 24（新增 CC 家族 motif）
    DomainCategory::Tir,     // 25（新增 TIR 家族 motif）
    DomainCategory::Tir,     // 26（新增 TIR 家族 motif）
    DomainCategory::Tir,     // 27（新增 TIR 家族 motif）
    DomainCategory::Tir,     // 28（新增 TIR 家族 motif）
];

/// RGB color table (index = motif id).
pub const RGB_COLORS: [[u8; 3]; 29] = [
    [0, 0, 0],
    [0, 255, 255],   // 1
    [0, 0, 255],     // 2
    [255, 0, 0],     // 3
    [255, 0, 255],   // 4
    [255, 255, 0],   // 5
    [0, 255, 0],     // 6
    [0, 128, 128],   // 7
    [68, 68, 68],    // 8
    [0, 128, 0],     // 9
    [192, 192, 192], // 10
    [128, 0, 128],   // 11
    [128, 128, 0],   // 12
    [0, 0, 128],     // 13
    [128, 0, 0],     // 14
    [255, 255, 255], // 15
    [0, 255, 255],   // 16
    [0, 0, 255],     // 17
    [255, 0, 0],     // 18
    [255, 0, 255],   // 19
    [255, 255, 0],   // 20
    [255, 128, 0],   // 21（RNL CC）
    [255, 192, 0],   // 22（RNL CC）
    [255, 64, 0],   // 23（RNL CC）
    [255, 96, 32],   // 24（RNL CC）
    [255, 160, 64],   // 25（RNL CC）
    [255, 128, 0],   // 26（RNL CC）
    [255, 192, 0],   // 27（RNL CC）
    [255, 64, 0],   // 28（RNL CC）
];

/// Enabled seed combinations.
pub const SEED_COMBINATIONS: &[&[MotifId]] = &[
    &[1, 6, 4],
    &[6, 4, 5],
    &[4, 5, 10],
    &[5, 10, 3],
    &[10, 3, 12],
    &[3, 12, 2],
    &[1, 4, 5],
    &[12, 2, 8],
    &[2, 8, 7],
    &[18, 15, 13],
    &[1, 6],
    // ---- 以下为 RNL/helper NLR 补丁新增：以 CC motif 开头、后接 NB-ARC motif ----
    &[18, 15],
    &[15, 13],
];

/// NLR signatures.
pub const SIGNATURES: &[&[MotifId]] = &[
    &[17, 16],
    &[1, 6],
    &[1, 6, 4],
    &[6, 4, 5],
    &[4, 5, 10],
    &[5, 10, 3],
    &[10, 3, 12],
    &[3, 12, 2],
    &[12, 2, 8],
    &[2, 8, 7],
    &[8, 7, 9],
    &[7, 9, 11],
    &[9, 11],
    &[11, 9],
    &[18, 15],
    &[15, 13],
    &[13, 1],
    &[1, 4, 5],
    // ---- RNL 补丁新增（用于 has_signature 预过滤）----
];

/// Consensus sequences for NB-ARC motifs.
/// Index 0 is a placeholder; only NB-ARC motifs (1,2,3,4,5,6,10,12) have values.
pub const CONSENSUS_SEQUENCES: [&str; 29] = [
    "",
    "PIWGMGGVGKTTLARAVYNDP",          // 1 (P-loop)
    "LKPCFLYCAIFPEDYMIDKNKLIWLWMAE",  // 2
    "CGGLPLAIKVWGGMLAGKQKT",          // 3
    "YLVVLDDVWDTDQWD",                // 4
    "NGSRIIITTRNKHVANYMCT",           // 5
    "HFDCRAWVCVSQQYDMKKVLRDIIQQVGG",  // 6
    "",                                // 7
    "",                                // 8
    "",                                // 9
    "LSHEESWQLFHQHAF",                // 10
    "",                                // 11
    "IMPVLRLSYHHLPYH",                // 12
    "", "", "", "", "", "", "", "",    // 13..20
    "",                                // 21（RNL CC，非 NB-ARC）
    "",                                // 22（RNL CC，非 NB-ARC）
    "",                                // 23（RNL CC，非 NB-ARC）
    "",                                // 24（RNL CC，非 NB-ARC）
    "",                                // 25（RNL CC，非 NB-ARC）
    "",                                // 26（RNL CC，非 NB-ARC）
    "",                                // 27（RNL CC，非 NB-ARC）
    "",                                // 28（RNL CC，非 NB-ARC）
];

/// Convert a motif id to a "motif_N" string.
#[inline]
pub fn motif_id_str(id: MotifId) -> String {
    format!("motif_{}", id)
}

/// NLR annotation signature definition (access interface for rank/category/color/seed/consensus/ploop).
///
/// 内置表都是编译期常量；额外类别（外部 motif 库用）放在 HashMap 里，因此本结构不再是零大小、
/// 也不再是 Copy——需要多处使用时 clone。
#[derive(Debug, Clone, Default)]
pub struct AnnotatorSignatureDefinition {
    /// 内置表之外的 motif 类别（外部库用命令行声明，例如 21=CC、25=TIR）。
    extra_categories: std::collections::HashMap<MotifId, DomainCategory>,
    /// 外部库额外声明的播种组合（形如 [21,4]）；含 >20 的 id 时用连续子串匹配。
    extra_seeds: Vec<Vec<MotifId>>,
    /// 外部库额外声明的 signature（预过滤用）。
    extra_signatures: Vec<Vec<MotifId>>,
}

/// Number of built-in motifs. Ids above this value are treated as library-specific motifs;
/// their categories default to `NA` and can be declared via `with_extra_categories()`.
pub const BUILTIN_MOTIF_COUNT: MotifId = 20;

impl AnnotatorSignatureDefinition {
    pub fn new() -> Self {
        AnnotatorSignatureDefinition {
            extra_categories: std::collections::HashMap::new(),
            extra_seeds: Vec::new(),
            extra_signatures: Vec::new(),
        }
    }

    /// 声明内置表之外 motif 的域类别（外部 mot.txt 带更多 motif 时使用）。
    pub fn with_extra_categories(mut self, m: std::collections::HashMap<MotifId, DomainCategory>) -> Self {
        self.extra_categories = m;
        self
    }

    /// 外部库额外声明的播种组合 / signature（内置表之外新 motif 的配套规则）。
    pub fn with_extra_rules(mut self, seeds: Vec<Vec<MotifId>>, signatures: Vec<Vec<MotifId>>) -> Self {
        self.extra_seeds = seeds;
        self.extra_signatures = signatures;
        self
    }

    /// Rank of the motif. 超出内置表时返回 15（最低档），不会 panic。
    #[inline]
    pub fn rank(&self, id: MotifId) -> u8 {
        RANKS.get(id as usize).copied().unwrap_or(15)
    }

    /// Domain category of the motif：内置表 → 命令行声明的额外类别 → NA。
    #[inline]
    pub fn category(&self, id: MotifId) -> DomainCategory {
        if let Some(c) = self.extra_categories.get(&id) {
            return *c;
        }
        CATEGORIES.get(id as usize).copied().unwrap_or(DomainCategory::Na)
    }

    /// Whether it is LRR.
    #[inline]
    pub fn is_lrr(&self, id: MotifId) -> bool {
        self.category(id) == DomainCategory::Lrr
    }

    /// Whether it is NBARC.
    #[inline]
    pub fn is_nbarc(&self, id: MotifId) -> bool {
        self.category(id) == DomainCategory::Nbarc
    }

    /// Whether it is P-loop.
    #[inline]
    pub fn is_ploop(&self, id: MotifId) -> bool {
        id == PLOOP_MOTIF
    }

    /// RGB color, returned as a "r,g,b" string.
    #[inline]
    pub fn color_rgb(&self, id: MotifId) -> String {
        let c = RGB_COLORS.get(id as usize).unwrap_or(&[128, 128, 128]);
        format!("{},{},{}", c[0], c[1], c[2])
    }

    /// Consensus sequence (empty string if none).
    #[inline]
    pub fn consensus(&self, id: MotifId) -> &'static str {
        CONSENSUS_SEQUENCES.get(id as usize).copied().unwrap_or("")
    }

    /// NB-ARC motif order sorted by rank ascending.
    /// Fixed as [1, 6, 4, 5, 10, 3, 12, 2].
    pub fn nbarc_motif_order(&self) -> Vec<MotifId> {
        let mut ids: Vec<MotifId> = (1..=BUILTIN_MOTIF_COUNT)
            .filter(|&i| self.is_nbarc(i))
            .collect();
        ids.sort_by_key(|&i| self.rank(i));
        ids
    }

    /// Whether the accumulated motif id sequence is a seed combination.
    ///
    /// 内置组合（只含 motif_1..motif_20）保持原版语义：**整串精确相等** `seq == combo`。
    /// 含"外部库新增 motif"（id > BUILTIN_MOTIF_COUNT）的组合改用**连续子串匹配**：
    /// 实测 helper NLR/RNL 区域的命中串形如 `[21, 23, 1, 6, 4, ...]`，若仍要求整串相等，
    /// 新 motif 永远无法播种（这也是原版 motif 法找不到 RNL 的原因之一）。
    pub fn is_seed(&self, seq: &[MotifId]) -> bool {
        let extra = self.extra_seeds.iter().map(|v| v.as_slice());
        SEED_COMBINATIONS
            .iter()
            .map(|v| *v)
            .chain(extra)
            .any(|s: &[MotifId]| {
            if s.iter().any(|&i| i > BUILTIN_MOTIF_COUNT) {
                s.len() <= seq.len() && seq.windows(s.len()).any(|w| w == s)
            } else {
                s == seq
            }
        })
    }

    /// 宽松播种（新增）：连续命中的 motif 串里 NB-ARC 类 motif 数 ≥ min_nbarc 即视为 seed。
    ///
    /// 动机：原版只认 11 个硬编码的精确 motif ID 组合（SEED_COMBINATIONS），而 helper NLR（RNL，
    /// 如 ADR1/NRG1）与部分分化 NLR 的 NB-ARC motif 命中顺序/组合不在表里（实测大豆 RNL 区域
    /// 命中 motif_2,motif_3,motif_11,motif_9 → 不出位点）。本模式按"域类别"而非"精确 ID 串"播种。
    pub fn is_seed_relaxed(&self, seq: &[MotifId], min_nbarc: usize) -> bool {
        let n_nbarc = seq
            .iter()
            .filter(|&&id| self.category(id) == DomainCategory::Nbarc)
            .count();
        n_nbarc >= min_nbarc
    }

    /// Whether the motif id sequence contains an NLR signature (a contiguous subsequence matches any signature).
    pub fn has_signature(&self, ids: &[MotifId]) -> bool {
        let extra = self.extra_signatures.iter().map(|v| v.as_slice());
        SIGNATURES
            .iter()
            .map(|v| *v)
            .chain(extra)
            .any(|sig: &[MotifId]| {
            sig.len() <= ids.len()
                && ids.windows(sig.len()).any(|w| w == sig)
        })
    }

    /// Derive the domain string from a motif list.
    /// Merges consecutive identical categories, skips NA and LINKER, outputs "NBARC-LRR" style.
    pub fn domain_string(&self, ids: &[MotifId]) -> String {
        let mut parts: Vec<&'static str> = Vec::new();
        let mut current: Option<DomainCategory> = None;
        for &id in ids {
            let cat = self.category(id);
            if cat == DomainCategory::Na || cat == DomainCategory::Linker {
                current = Some(cat); // skipped but updates continuity.
                continue;
            }
            if current != Some(cat) {
                parts.push(cat.as_str());
                current = Some(cat);
            }
        }
        parts.join("-")
    }
}

/// Signature pre-filter definition. The actual logic is merged into `AnnotatorSignatureDefinition::has_signature`.
#[derive(Debug, Clone, Copy, Default)]
pub struct SignatureDefinition;

impl SignatureDefinition {
    pub fn new() -> Self {
        SignatureDefinition
    }

    /// Whether any signature matches (delegates to the signature table).
    #[inline]
    pub fn has_signature(&self, ids: &[MotifId]) -> bool {
        AnnotatorSignatureDefinition::new().has_signature(ids)
    }
}
