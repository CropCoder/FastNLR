//! nlr-assemble integration tests: validate the three-step assembly.

use nlr_assemble::{assemble, AssembleParams};
use nlr_core::motif::Motif;
use nlr_core::signature_def::AnnotatorSignatureDefinition;
use nlr_core::strand::Strand;

fn def() -> AnnotatorSignatureDefinition {
    AnnotatorSignatureDefinition::new()
}

/// Build a motif with DNA coordinates already set.
fn dna_motif(id: u8, start: u64, strand: Strand) -> Motif {
    let mut m = Motif::new_protein(id, "chr1".to_string(), 1, "AAA".to_string(), 1e-6);
    m.set_dna("chr1".to_string(), start, 20000, 0, strand);
    m
}

#[test]
fn empty_motifs_no_nlrs() {
    let nlrs = assemble("chr1", vec![], &AssembleParams::default(), &def());
    assert!(nlrs.is_empty());
}

#[test]
fn seed_combination_produces_nlr() {
    // Build a seed combination 1,6,4 (within 100bp), should produce one NLR.
    let motifs = vec![
        dna_motif(1, 1000, Strand::Forward),
        dna_motif(6, 1100, Strand::Forward),
        dna_motif(4, 1200, Strand::Forward),
    ];
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &def());
    assert_eq!(nlrs.len(), 1);
    assert_eq!(nlrs[0].name, "chr1_nlr1");
}

#[test]
fn no_seed_no_nlr() {
    // Motifs without a seed combination produce no NLR.
    let motifs = vec![
        dna_motif(1, 1000, Strand::Forward),
        dna_motif(9, 20000, Strand::Forward), // LRR, but no NBARC seed
    ];
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &def());
    assert!(nlrs.is_empty());
}

#[test]
fn overlapping_seeds_do_not_duplicate_motifs() {
    // 间距 300/300/150/150 bp（均 <= distance_within_motif_combination = 500）会形成互相重叠的种子：
    // [1,6] 命中后，[6,4,5]、[4,5,10] 又从各自起点重复并入相同的 motif。
    let motifs = vec![
        dna_motif(1, 0, Strand::Forward),
        dna_motif(6, 300, Strand::Forward),
        dna_motif(4, 600, Strand::Forward),
        dna_motif(5, 750, Strand::Forward),
        dna_motif(10, 900, Strand::Forward),
    ];
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &def());
    assert_eq!(nlrs.len(), 1);

    let ids: Vec<u8> = nlrs[0].motifs.iter().map(|m| m.id).collect();
    assert_eq!(ids, vec![1, 6, 4, 5, 10], "位点内不应出现重复 motif");

    let mut keys: Vec<(u8, u64)> = nlrs[0].motifs.iter().map(|m| (m.id, m.dna_start)).collect();
    let total = keys.len();
    keys.sort();
    keys.dedup();
    assert_eq!(keys.len(), total, "同一 (id, dna_start) 只能出现一次");
}
