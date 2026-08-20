//! NLR loci (a collection of motifs on one sequence) type.
//!
//! Replicates Java `MotifList`:
//! - Sorting (by Motif's Ord; forward strand ascending / reverse strand descending);
//! - `is_complete_nlr` (contains P-loop and LRR);
//! - `can_be_merged_with` (same strand / nearest of the 4 dna_start combinations within threshold / rank monotonic);
//! - `remove_redundant_motifs` (dedup by id+start+end).

use crate::motif::Motif;
use crate::signature_def::{AnnotatorSignatureDefinition, MotifId};
use crate::strand::Strand;

/// An NLR loci (an ordered collection of motifs on one sequence).
#[derive(Debug, Clone, PartialEq)]
pub struct MotifList {
    /// Loci name (e.g. `chr1_nlr1`).
    pub name: String,
    /// Ordered motif collection.
    pub motifs: Vec<Motif>,
}

impl MotifList {
    pub fn new(name: String, motifs: Vec<Motif>) -> Self {
        let mut list = MotifList { name, motifs };
        list.sort();
        list
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

    /// Whether it is a complete NLR: contains P-loop and LRR (Java `isCompleteNLR`).
    pub fn is_complete_nlr(&self, def: &AnnotatorSignatureDefinition) -> bool {
        let has_ploop = self.motifs.iter().any(|m| def.is_ploop(m.id));
        let has_lrr = self.motifs.iter().any(|m| def.is_lrr(m.id));
        has_ploop && has_lrr
    }

    /// Whether it contains a stop codon (Java `hasStopCodon`).
    pub fn has_stop_codon(&self) -> bool {
        self.motifs.iter().any(|m| m.has_stop())
    }

    /// Whether it can be merged with another loci (Java `canBeMergedWith`).
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
        // Replicates Java: HashSet<Motif> (equals = id + dnaStart) dedup then sort;
        // during iteration a motif passes only if `rank > lastRank || isLRR`, otherwise
        // (rank <= lastRank and not LRR) returns false.
        let mut merged: Vec<&Motif> = self.motifs.iter().chain(other.motifs.iter()).collect();
        // Dedup by (id, dna_start) (equivalent to Java Motif.equals).
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

    /// Whether it contains a used motif (Java `containsUsedMotif`, used for assembly dedup).
    ///
    /// Note: assembly dedup is actually implemented in nlr-assemble keyed by `(id, dna_start)`
    /// (equivalent to Java `Motif.equals`); this method is retained only for API completeness
    /// and is not part of the current assembly flow.
    pub fn contains_used_motif(&self, used: &std::collections::HashSet<(u8, u64)>) -> bool {
        self.motifs
            .iter()
            .any(|m| used.contains(&(m.id, m.dna_start)))
    }

    /// Dedup by id+start+end (Java `removeRedundantMotifs`, keeps first occurrence).
    pub fn remove_redundant_motifs(&mut self) {
        let mut seen = std::collections::HashSet::new();
        self.motifs.retain(|m| {
            let key = (m.id, m.dna_start, m.dna_end);
            seen.insert(key)
        });
    }

    /// Motif id list (comma-joined, Java `getMotifListString`).
    pub fn motif_list_string(&self) -> String {
        self.motifs
            .iter()
            .map(|m| crate::signature_def::motif_id_str(m.id))
            .collect::<Vec<_>>()
            .join(",")
    }

    /// Domain string (Java `getDomainString`).
    pub fn domain_string(&self, def: &AnnotatorSignatureDefinition) -> String {
        let ids: Vec<MotifId> = self.motifs.iter().map(|m| m.id).collect();
        def.domain_string(&ids)
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
    /// Replicates Java `MotifList.compareTo`: same strand forward ascending / reverse descending; on equal start, the longer one comes first.
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
        match start_cmp {
            Ordering::Equal => other.motifs.len().cmp(&self.motifs.len()),
            o => o,
        }
    }
}
