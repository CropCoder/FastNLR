//! nlr-core unit tests: verifies coordinate mapping, sorting, seed/signature tables, merging and completeness.

#[cfg(test)]
mod tests {
    use nlr_core::motif::Motif;
    use nlr_core::motif_list::MotifList;
    use nlr_core::signature_def::AnnotatorSignatureDefinition;
    use nlr_core::strand::Strand;

    fn def() -> AnnotatorSignatureDefinition {
        AnnotatorSignatureDefinition::new()
    }

    fn protein_motif(id: u8, position: u64, seq: &str, pvalue: f64) -> Motif {
        Motif::new_protein(id, "seq".to_string(), position, seq.to_string(), pvalue)
    }

    #[test]
    fn set_dna_forward_strand() {
        let mut m = protein_motif(1, 1, "PIWG", 1e-6);
        // Sequence length 4; forward strand frame 0 offset 100.
        m.set_dna("chr1".into(), 100, 20000, 0, Strand::Forward);
        // dna_start = (1-1)*3 + 0 + 100 = 100; dna_end = (1+4-1)*3 + 0 + 100 = 112.
        assert_eq!(m.dna_start, 100);
        assert_eq!(m.dna_end, 112);
        assert!(m.dna_parameters_set);
    }

    #[test]
    fn set_dna_reverse_strand() {
        let mut m = protein_motif(1, 1, "PIWG", 1e-6);
        m.set_dna("chr1".into(), 100, 20000, 0, Strand::Reverse);
        // dna_start = 100 + 20000 - ((1+4-1)*3 + 0) = 20100 - 12 = 20088;
        // dna_end   = 100 + 20000 - ((1-1)*3 + 0) = 20100.
        assert_eq!(m.dna_start, 20088);
        assert_eq!(m.dna_end, 20100);
        assert!(m.dna_start < m.dna_end); // reverse strand coordinates are not inverted; start is always leftmost
    }

    #[test]
    fn seed_and_signature_tables() {
        let d = def();
        assert!(d.is_seed(&[1, 6, 4]));
        assert!(d.is_seed(&[1, 6]));
        assert!(!d.is_seed(&[6, 1])); // order-sensitive
        assert!(d.has_signature(&[1, 6, 4, 5]));
        assert!(!d.has_signature(&[1, 2, 3]));
    }

    #[test]
    fn nbarc_order_matches_java() {
        let d = def();
        assert_eq!(d.nbarc_motif_order(), vec![1, 6, 4, 5, 10, 3, 12, 2]);
    }

    #[test]
    fn domain_string_skips_linker_na() {
        let d = def();
        // TIR(13) NBARC(1) LINKER(7) LRR(9) -> "TIR-NBARC-LRR"
        let s = d.domain_string(&[13, 1, 7, 9]);
        assert_eq!(s, "TIR-NBARC-LRR");
    }

    #[test]
    fn is_complete_nlr() {
        let d = def();
        let complete = MotifList::new(
            "x".into(),
            vec![
                protein_motif(1, 1, "P", 1e-6), // P-loop
                protein_motif(9, 20, "L", 1e-6), // LRR
            ],
        );
        assert!(complete.is_complete_nlr(&d));

        let incomplete = MotifList::new("x".into(), vec![protein_motif(1, 1, "P", 1e-6)]);
        assert!(!incomplete.is_complete_nlr(&d));
    }

    #[test]
    fn can_be_merged_with_rank_consistency() {
        let d = def();
        // Two same-strand loci, rank monotonic (1 < 6), and within distance.
        let a = MotifList::new("a".into(), vec![protein_motif(1, 1, "A", 1e-6)]);
        let b = MotifList::new("b".into(), vec![protein_motif(6, 2, "B", 1e-6)]);
        // Set DNA so they are on the same strand and close in distance.
        let mut am = a.motifs[0].clone();
        let mut bm = b.motifs[0].clone();
        am.set_dna("chr".into(), 0, 20000, 0, Strand::Forward);
        bm.set_dna("chr".into(), 100, 20000, 0, Strand::Forward);
        let a2 = MotifList::new("a".into(), vec![am]);
        let b2 = MotifList::new("b".into(), vec![bm]);
        assert!(a2.can_be_merged_with(&b2, 500, &d));
    }
}
