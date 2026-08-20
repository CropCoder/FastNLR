//! Regression tests: verifies the 6 algorithm differences fixed in this round.

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;
use nlr_core::strand::Strand;

fn def() -> AnnotatorSignatureDefinition {
    AnnotatorSignatureDefinition::new()
}

fn dna_motif(id: u8, dna_seq_id: &str, protein_seq_id: &str, start: u64, strand: Strand) -> Motif {
    let mut m = Motif::new_protein(id, protein_seq_id.to_string(), 1, "AAA".to_string(), 1e-6);
    m.set_dna(dna_seq_id.to_string(), start, 20000, 0, strand);
    m
}

#[test]
fn motif_cmp_groups_by_dna_sequence_id() {
    // Key bug: previously protein_sequence_id (fragment id) was used to determine the same group,
    // so same-chromosome motifs were never sorted by coordinate.
    // Now dna_sequence_id is used to determine the same group; within a group, sort by dna_start.
    let m1 = dna_motif(1, "chr1", "chr1_50000_frame+0", 1000, Strand::Forward);
    let m2 = dna_motif(6, "chr1", "chr1_0_frame+0", 500, Strand::Forward);
    // Different protein_sequence_id (different fragments) but same dna_sequence_id (same chromosome).
    // Correct order: m2(500) < m1(1000), because sorted by dna_start.
    let mut v = vec![m1.clone(), m2.clone()];
    v.sort();
    assert_eq!(v[0].dna_start, 500, "should sort by dna_start, got {:?}", v.iter().map(|m| m.dna_start).collect::<Vec<_>>());
}

#[test]
fn can_be_merged_rank_equal_non_lrr_rejected() {
    // Key bug: equal rank and non-LRR should reject (previously missed).
    let d = def();
    let mut m1 = dna_motif(1, "chr1", "p1", 100, Strand::Forward); // rank 4
    let mut m2 = dna_motif(6, "chr1", "p2", 200, Strand::Forward); // rank 5
    let mut m3 = dna_motif(1, "chr1", "p3", 300, Strand::Forward); // rank 4 (same rank as m1, non-LRR)
    m1.set_dna("chr1".into(), 100, 20000, 0, Strand::Forward);
    m2.set_dna("chr1".into(), 200, 20000, 0, Strand::Forward);
    m3.set_dna("chr1".into(), 300, 20000, 0, Strand::Forward);

    let a = MotifList::new("a".into(), vec![m1, m2]); // rank 4,5 monotonic
    let b = MotifList::new("b".into(), vec![m3]); // rank 4
    // After merge rank 4,5,4 -> 4 then 4 again (non-LRR), should reject.
    assert!(!a.can_be_merged_with(&b, 10000, &d));
}

#[test]
fn span_reverse_strand() {
    // Reverse strand motif dna_start is always leftmost; span should be min/max.
    let mut m = Motif::new_protein(1, "p".into(), 1, "AAA".into(), 1e-6);
    m.set_dna("chr1".into(), 100, 20000, 0, Strand::Reverse);
    // Reverse strand: dna_start = 100 + 20000 - ((1+3-1)*3) = 20100 - 9 = 20091; dna_end = 100 + 20000 - 0 = 20100
    assert!(m.dna_start < m.dna_end, "reverse strand dna_start should be < dna_end");
    let list = MotifList::new("x".into(), vec![m]);
    let (start, end) = list.span();
    assert!(start < end);
}

#[test]
fn elongate_uses_dynamic_first_motif() {
    // Verify that add_motif triggers a re-sort and first_motif/last_motif return dynamic results.
    let mut m1 = dna_motif(6, "chr1", "p1", 500, Strand::Forward); // rank 5
    let mut m4 = dna_motif(4, "chr1", "p2", 300, Strand::Forward); // rank 6
    m1.set_dna("chr1".into(), 500, 20000, 0, Strand::Forward);
    m4.set_dna("chr1".into(), 300, 20000, 0, Strand::Forward);

    let mut list = MotifList::new("seed".into(), vec![m1]);
    assert_eq!(list.first_motif().dna_start, 500);
    // After adding a motif with dna_start=300, first_motif should become 300.
    list.add_motif(m4);
    assert_eq!(list.first_motif().dna_start, 300, "after add_motif, first_motif should update dynamically");
}
