//! NLR loci (a collection of motifs on one sequence) type.
//!
//! Key behavior:
//! - Sorting (by Motif's Ord; forward strand ascending / reverse strand descending);
//! - `is_complete_nlr` (contains P-loop and LRR);
//! - `can_be_merged_with` (same strand / nearest of the 4 dna_start combinations within threshold / rank monotonic);
//! - `remove_redundant_motifs` (dedup by id+start+end).

use crate::motif::Motif;
use crate::signature_def::{AnnotatorSignatureDefinition, DomainCategory, DNA_REPEAT_MIN_BLOCK};
use crate::strand::Strand;

/// An NLR loci (an ordered collection of motifs on one sequence).
#[derive(Debug, Clone, PartialEq)]
pub struct MotifList {
    /// Loci name (e.g. `chr1_nlr1`).
    pub name: String,
    /// Ordered motif collection.
    pub motifs: Vec<Motif>,
    /// HMM 证据补充的 N 端结构域（仅 DNA 模式会设置）。
    ///
    /// DNA 模式的 N 端域判定依赖 PWM motif，而内置库对 CC 的覆盖有限（尤其没有 RPW8 型 CC），
    /// 会把 RNL/helper NLR 与部分 CNL 误判为 `NBARC-LRR`。因此在组装完成后对位点 ORF 再跑一次
    /// 内置 HMM，若发现 PWM 未给出的 CC/TIR，就记在这里，由 [`MotifList::domain_string`] 放到
    /// 域串最前面。`None` 表示沿用 PWM 判定。
    pub hmm_nterm: Option<DomainCategory>,
    /// HMM 证据给出的生物学家族（TNL/RNL/CNL/NB-only/TIR-only）。
    pub hmm_family: Option<String>,
}

impl MotifList {
    pub fn new(name: String, motifs: Vec<Motif>) -> Self {
        let mut list = MotifList {
            name,
            motifs,
            hmm_nterm: None,
            hmm_family: None,
        };
        list.sort();
        list
    }

    /// 由 HMM 证据补充 N 端结构域。
    pub fn set_hmm_nterm(&mut self, category: DomainCategory) {
        self.hmm_nterm = Some(category);
    }

    /// Store the family inferred from HMM evidence for DNA-mode reporting.
    pub fn set_hmm_family(&mut self, family: impl Into<String>) {
        self.hmm_family = Some(family.into());
    }

    /// Sort (forward strand ascending, reverse strand descending).
    pub fn sort(&mut self) {
        self.motifs.sort();
    }

    /// First motif.
    #[inline]
    pub fn first_motif(&self) -> &Motif {
        &self.motifs[0]
    }

    /// Last motif.
    #[inline]
    pub fn last_motif(&self) -> &Motif {
        &self.motifs[self.motifs.len() - 1]
    }

    /// Strand (determined by the first motif).
    #[inline]
    pub fn strand(&self) -> Strand {
        self.first_motif().strand
    }

    /// Whether on the forward strand.
    #[inline]
    pub fn is_forward(&self) -> bool {
        self.strand().is_forward()
    }

    /// Sequence name (DNA sequence id, determined by the first motif).
    #[inline]
    pub fn sequence_name(&self) -> &str {
        self.first_motif()
            .dna_sequence_id
            .as_deref()
            .unwrap_or(&self.first_motif().protein_sequence_id)
    }

    /// Append a motif and re-sort.
    pub fn add_motif(&mut self, motif: Motif) {
        self.motifs.push(motif);
        self.sort();
    }

    /// Append a set of motifs and re-sort.
    pub fn add_motifs(&mut self, motifs: &[Motif]) {
        self.motifs.extend(motifs.iter().cloned());
        self.sort();
    }

    /// Whether it is a complete NLR: contains P-loop and LRR.
    pub fn is_complete_nlr(&self, def: &AnnotatorSignatureDefinition) -> bool {
        let has_ploop = self.motifs.iter().any(|m| def.is_ploop(m.id));
        let has_lrr = self.motifs.iter().any(|m| def.is_lrr(m.id));
        has_ploop && has_lrr
    }

    /// Whether it contains a stop codon.
    pub fn has_stop_codon(&self) -> bool {
        self.motifs.iter().any(|m| m.has_stop())
    }

    /// Whether it can be merged with another loci.
    ///
    /// Three conditions:
    /// 1. Same strand;
    /// 2. At least one of the 4 dna_start endpoint combinations is within the threshold;
    /// 3. After merging, rank is monotonically non-decreasing along the translation direction (LRR may tie).
    pub fn can_be_merged_with(
        &self,
        other: &MotifList,
        distance: u64,
        def: &AnnotatorSignatureDefinition,
    ) -> bool {
        if self.strand() != other.strand() {
            return false;
        }

        // Distance: 4 dna_start combinations.
        let a_first = self.first_motif().dna_start;
        let a_last = self.last_motif().dna_start;
        let b_first = other.first_motif().dna_start;
        let b_last = other.last_motif().dna_start;
        let all_far = (a_first.abs_diff(b_first) > distance)
            && (a_first.abs_diff(b_last) > distance)
            && (a_last.abs_diff(b_first) > distance)
            && (a_last.abs_diff(b_last) > distance);
        if all_far {
            return false;
        }

        // Rank consistency: after merge, dedup, and sort, rank is monotonically non-decreasing (LRR exception).
        // Deduplicate by motif identity and DNA start, then sort.
        // during iteration a motif passes only if `rank > lastRank || isLRR`, otherwise
        // (rank <= lastRank and not LRR) returns false.
        let mut merged: Vec<&Motif> = self.motifs.iter().chain(other.motifs.iter()).collect();
        // Deduplicate by (id, dna_start).
        let mut seen = std::collections::HashSet::new();
        merged.retain(|m| seen.insert((m.id, m.dna_start)));
        merged.sort();
        if merged.is_empty() {
            return true;
        }
        let mut last_rank = def.rank(merged[0].id);
        for m in merged.iter().skip(1) {
            let rank = def.rank(m.id);
            if rank > last_rank || def.is_lrr(m.id) {
                last_rank = rank;
            } else {
                return false;
            }
        }
        true
    }

    /// Whether it contains a used motif (used for assembly dedup).
    ///
    /// Note: assembly dedup is actually implemented in nlr-assemble keyed by `(id, dna_start)`
    /// This method is retained only for API completeness and is not part of the current assembly flow.
    pub fn contains_used_motif(&self, used: &std::collections::HashSet<(u8, u64)>) -> bool {
        self.motifs
            .iter()
            .any(|m| used.contains(&(m.id, m.dna_start)))
    }

    /// Dedup by id+start+end, keeping the first occurrence.
    pub fn remove_redundant_motifs(&mut self) {
        let mut seen = std::collections::HashSet::new();
        self.motifs.retain(|m| {
            let key = (m.id, m.dna_start, m.dna_end);
            seen.insert(key)
        });
    }

    /// Motif id list (comma-joined).
    /// Domain string.
    pub fn domain_string(&self, def: &AnnotatorSignatureDefinition) -> String {
        let mut cats: Vec<DomainCategory> =
            self.motifs.iter().map(|m| def.category(m.id)).collect();
        // HMM 补充的 N 端域位于域串最前（仅当 PWM 通路未能给出 N 端域时才会被设置）。
        if let Some(c) = self.hmm_nterm {
            cats.insert(0, c);
        }
        crate::signature_def::canonical_domain_string(&cats, DNA_REPEAT_MIN_BLOCK).1
    }

    /// Biological family, separating RNL from CNL when they share `CC-NBARC-LRR`.
    pub fn family<'a>(&'a self, def: &AnnotatorSignatureDefinition) -> &'a str {
        if let Some(family) = self.hmm_family.as_deref() {
            return family;
        }
        let class = self.domain_string(def);
        if class.starts_with("TIR") {
            if class.contains("NBARC") {
                "TNL"
            } else {
                "TIR-only"
            }
        } else if class.starts_with("CC") {
            if self.motifs.iter().any(|m| matches!(m.id, 21 | 23)) {
                "RNL"
            } else {
                "CNL"
            }
        } else if class.contains("NBARC") {
            "NB-only"
        } else {
            "Other"
        }
    }

    /// Loci span (0-based, min/max of the four first/last values).
    pub fn span(&self) -> (u64, u64) {
        let a = self.first_motif();
        let b = self.last_motif();
        let start = a.dna_start.min(a.dna_end).min(b.dna_start).min(b.dna_end);
        let end = a.dna_start.max(a.dna_end).max(b.dna_start).max(b.dna_end);
        (start, end)
    }
}

impl PartialOrd for MotifList {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for MotifList {}

impl Ord for MotifList {
    /// Sort order: forward strand ascending / reverse strand descending; on equal start, the longer loci comes first.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        use std::cmp::Ordering;
        let s = self.first_motif();
        let o = other.first_motif();
        match (s.strand, o.strand) {
            (Strand::Forward, Strand::Reverse) => return Ordering::Less,
            (Strand::Reverse, Strand::Forward) => return Ordering::Greater,
            _ => {}
        }
        let start_cmp = match s.strand {
            Strand::Forward => s.dna_start.cmp(&o.dna_start),
            Strand::Reverse => o.dna_start.cmp(&s.dna_start),
        };
        if start_cmp != Ordering::Equal {
            return start_cmp;
        }
        let self_key: Vec<_> = self
            .motifs
            .iter()
            .map(|m| (m.id, m.dna_start, m.dna_end, m.frame, m.position))
            .collect();
        let other_key: Vec<_> = other
            .motifs
            .iter()
            .map(|m| (m.id, m.dna_start, m.dna_end, m.frame, m.position))
            .collect();
        other
            .motifs
            .len()
            .cmp(&self.motifs.len())
            .then_with(|| self_key.cmp(&other_key))
            .then_with(|| {
                self.hmm_nterm
                    .map(DomainCategory::as_str)
                    .cmp(&other.hmm_nterm.map(DomainCategory::as_str))
            })
            .then_with(|| self.hmm_family.cmp(&other.hmm_family))
            .then_with(|| self.name.cmp(&other.name))
    }
}
