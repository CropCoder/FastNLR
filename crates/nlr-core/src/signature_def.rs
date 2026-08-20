//! NLR static rule tables and signature definitions.
//!
//! Faithfully replicates the hard-coded rules of Java `AnnotatorSignatureDefinition` and `SignatureDefinition`:
//! - rank, category (domain category), and RGB color for each of the 20 motifs;
//! - 11 seed combinations (findSeeds targets);
//! - 18 "NLR signatures" (contiguous motif sequence patterns used for pre-filtering);
//! - consensus sequences for 8 NB-ARC motifs;
//! - P-loop motif id (motif_1).

/// Motif id: 1..=20 (corresponding to Java's "motif_1" .. "motif_20").
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
    /// Construct from a Java category string.
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

    /// Convert back to a Java category string (for output).
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

/// P-loop motif id (Java `ploopMotif = "motif_1"`).
pub const PLOOP_MOTIF: MotifId = 1;

/// Rank table (index = motif id; `RANKS[0]` is an unused placeholder).
/// Corresponds to Java `loadDefaultMotifRanks`.
pub const RANKS: [u8; 21] = [
    0, 4, 11, 9, 6, 7, 5, 13, 12, 14, 8, 14, 10, 3, 3, 2, 2, 1, 1, 14, 15,
];

/// Category table (index = motif id).
/// Corresponds to Java `loadDefaultMotifCategories`.
pub const CATEGORIES: [DomainCategory; 21] = [
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
];

/// RGB color table (index = motif id).
/// Corresponds to Java `loadDefaultMotifRgbColors`.
pub const RGB_COLORS: [[u8; 3]; 21] = [
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
];

/// 11 seed combinations (enabled; the uncommented subset of Java `loadDefaultMotifIDCombinations`).
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
];

/// 18 "NLR signatures" (Java `SignatureDefinition.loadDefaultSignature`).
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
];

/// Consensus sequences for 8 NB-ARC motifs (Java `loadDefaultDefaultMotifSequences`).
/// Index 0 is a placeholder; only NB-ARC motifs (1,2,3,4,5,6,10,12) have values.
pub const CONSENSUS_SEQUENCES: [&str; 21] = [
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
];

/// Convert a motif id to a "motif_N" string.
#[inline]
pub fn motif_id_str(id: MotifId) -> String {
    format!("motif_{}", id)
}

/// NLR annotation signature definition (access interface for rank/category/color/seed/consensus/ploop).
///
/// All data is compile-time const; the struct is zero-sized (unit struct).
#[derive(Debug, Clone, Copy, Default)]
pub struct AnnotatorSignatureDefinition;

impl AnnotatorSignatureDefinition {
    pub fn new() -> Self {
        AnnotatorSignatureDefinition
    }

    /// Rank of the motif.
    #[inline]
    pub fn rank(&self, id: MotifId) -> u8 {
        RANKS[id as usize]
    }

    /// Domain category of the motif.
    #[inline]
    pub fn category(&self, id: MotifId) -> DomainCategory {
        CATEGORIES[id as usize]
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
        let c = &RGB_COLORS[id as usize];
        format!("{},{},{}", c[0], c[1], c[2])
    }

    /// Consensus sequence (empty string if none).
    #[inline]
    pub fn consensus(&self, id: MotifId) -> &'static str {
        CONSENSUS_SEQUENCES[id as usize]
    }

    /// NB-ARC motif order sorted by rank ascending (Java `getNbarcMotifOrder`).
    /// Fixed as [1, 6, 4, 5, 10, 3, 12, 2].
    pub fn nbarc_motif_order(&self) -> Vec<MotifId> {
        let mut ids: Vec<MotifId> = (1..=20)
            .filter(|&i| self.is_nbarc(i))
            .collect();
        ids.sort_by_key(|&i| self.rank(i));
        ids
    }

    /// Whether the accumulated motif id sequence is a seed combination (exact match).
    pub fn is_seed(&self, seq: &[MotifId]) -> bool {
        SEED_COMBINATIONS.iter().any(|s| *s == seq)
    }

    /// Whether the motif id sequence contains an NLR signature (a contiguous subsequence matches any signature).
    pub fn has_signature(&self, ids: &[MotifId]) -> bool {
        SIGNATURES.iter().any(|sig| {
            sig.len() <= ids.len()
                && ids.windows(sig.len()).any(|w| w == *sig)
        })
    }

    /// Derive the domain string from a motif list (Java `getDomainString`).
    /// Merges consecutive identical categories, skips NA and LINKER, outputs "NBARC-LRR" style.
    pub fn domain_string(&self, ids: &[MotifId]) -> String {
        let mut parts: Vec<&'static str> = Vec::new();
        let mut current: Option<DomainCategory> = None;
        for &id in ids {
            let cat = self.category(id);
            if cat == DomainCategory::Na || cat == DomainCategory::Linker {
                current = Some(cat); // skipped but updates continuity (equivalent to Java's unconditional currentCategory update)
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

/// Signature pre-filter definition (equivalent to Java `SignatureDefinition`; actual logic has been merged into `AnnotatorSignatureDefinition::has_signature`).
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
