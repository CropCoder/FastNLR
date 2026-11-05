//! Locus-level structural filter for high-precision DNA-mode output.
//!
//! The filter runs after HMM N-terminal supplementation and before any locus
//! output is written. It keeps loci with a complete LRR architecture and
//! applies stricter evidence thresholds to partial loci (TIR-only or NB-ARC
//! fragments without LRR). It only uses evidence inside the assembled locus,
//! never an external gold standard.

use std::collections::HashSet;

use crate::motif_list::MotifList;
use crate::signature_def::{AnnotatorSignatureDefinition, DomainCategory};

/// Parameters for the high-precision locus filter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocusFilterConfig {
    /// Keep a TIR-only locus only when it has at least this many distinct TIR motifs.
    pub min_tir_motifs: usize,
    /// Keep a TIR-only locus only when its span is at least this many bp.
    pub min_tir_span: u64,
    /// When true, an NB-ARC partial locus (no LRR) must contain the P-loop motif (motif_1).
    pub require_ploop: bool,
    /// Keep an NB-ARC partial locus only when it has at least this many distinct NB-ARC motifs.
    pub min_nbarc_motifs: usize,
    /// Keep an NB-ARC partial locus only when its span is at least this many bp.
    pub min_nbarc_span: u64,
    /// Global minimum locus span (0 disables this check).
    pub min_locus_span: u64,
}

impl LocusFilterConfig {
    /// Default high-precision configuration.
    ///
    /// These values are selected by the three-dataset regression grid and are
    /// intended to reduce false positives without losing any gold-standard gene.
    pub const fn high_precision_default() -> Self {
        LocusFilterConfig {
            min_tir_motifs: 2,
            min_tir_span: 150,
            require_ploop: false,
            min_nbarc_motifs: 0,
            min_nbarc_span: 0,
            min_locus_span: 0,
        }
    }

    /// Disabled filter: keep every assembled locus (legacy behavior).
    pub const fn disabled() -> Self {
        LocusFilterConfig {
            min_tir_motifs: 0,
            min_tir_span: 0,
            require_ploop: false,
            min_nbarc_motifs: 0,
            min_nbarc_span: 0,
            min_locus_span: 0,
        }
    }
}

impl Default for LocusFilterConfig {
    fn default() -> Self {
        Self::high_precision_default()
    }
}

/// Keep or drop one locus according to the high-precision rules.
pub fn should_keep_locus(
    list: &MotifList,
    def: &AnnotatorSignatureDefinition,
    cfg: &LocusFilterConfig,
) -> bool {
    let (start, end) = list.span();
    let span = end.saturating_sub(start).saturating_add(1);

    if cfg.min_locus_span > 0 && span < cfg.min_locus_span {
        return false;
    }

    let mut distinct_tir = HashSet::new();
    let mut distinct_nbarc = HashSet::new();
    let mut has_lrr = false;
    let mut has_ploop = false;
    let mut has_tir = false;
    let mut has_nbarc = false;

    for motif in &list.motifs {
        if def.is_lrr(motif.id) {
            has_lrr = true;
        }
        if def.is_ploop(motif.id) {
            has_ploop = true;
        }
        match def.category(motif.id) {
            DomainCategory::Tir => {
                has_tir = true;
                distinct_tir.insert(motif.id);
            }
            DomainCategory::Nbarc => {
                has_nbarc = true;
                distinct_nbarc.insert(motif.id);
            }
            _ => {}
        }
    }

    if has_lrr {
        return true;
    }

    if has_nbarc {
        let ploop_ok = !cfg.require_ploop || has_ploop;
        return ploop_ok
            && distinct_nbarc.len() >= cfg.min_nbarc_motifs
            && span >= cfg.min_nbarc_span;
    }

    if has_tir {
        return distinct_tir.len() >= cfg.min_tir_motifs && span >= cfg.min_tir_span;
    }

    false
}

/// Filter a locus list in place.
pub fn filter_loci(
    nlrs: &mut Vec<MotifList>,
    def: &AnnotatorSignatureDefinition,
    cfg: &LocusFilterConfig,
) {
    nlrs.retain(|list| should_keep_locus(list, def, cfg));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::motif::Motif;
    use crate::strand::Strand;

    fn motif(id: u8, dna_start: u64) -> Motif {
        Motif {
            id,
            protein_sequence_id: "seq1".to_string(),
            position: 1,
            protein_sequence: "A".to_string(),
            pvalue: 0.0,
            score: 0.0,
            dna_sequence_id: Some("seq1".to_string()),
            dna_start,
            dna_end: dna_start,
            strand: Strand::Forward,
            frame: 0,
            dna_parameters_set: true,
        }
    }

    fn def() -> AnnotatorSignatureDefinition {
        AnnotatorSignatureDefinition::new()
    }

    fn strict() -> LocusFilterConfig {
        LocusFilterConfig {
            min_tir_motifs: 3,
            min_tir_span: 300,
            require_ploop: true,
            min_nbarc_motifs: 4,
            min_nbarc_span: 500,
            min_locus_span: 0,
        }
    }

    #[test]
    fn keeps_lrr_locus() {
        let list = MotifList::new(
            "s".into(),
            vec![motif(1, 100), motif(9, 200), motif(11, 300)],
        );
        assert!(should_keep_locus(
            &list,
            &def(),
            &LocusFilterConfig::default()
        ));
    }

    #[test]
    fn drops_short_tir_only_locus() {
        let list = MotifList::new("s".into(), vec![motif(13, 100), motif(15, 120)]);
        assert!(!should_keep_locus(
            &list,
            &def(),
            &LocusFilterConfig::default()
        ));
    }

    #[test]
    fn drops_nbarc_partial_without_ploop() {
        let list = MotifList::new(
            "s".into(),
            vec![motif(2, 100), motif(3, 130), motif(4, 160), motif(5, 190)],
        );
        assert!(!should_keep_locus(&list, &def(), &strict()));
    }

    #[test]
    fn keeps_nbarc_partial_with_ploop_and_enough_motifs() {
        let list = MotifList::new(
            "s".into(),
            vec![
                motif(1, 100),
                motif(2, 200),
                motif(3, 300),
                motif(4, 400),
                motif(5, 700),
            ],
        );
        assert!(should_keep_locus(&list, &def(), &strict()));
    }

    #[test]
    fn keeps_nbarc_partial_under_default_filter() {
        let list = MotifList::new(
            "s".into(),
            vec![motif(2, 100), motif(3, 130), motif(4, 160)],
        );
        assert!(should_keep_locus(
            &list,
            &def(),
            &LocusFilterConfig::default()
        ));
    }
}
