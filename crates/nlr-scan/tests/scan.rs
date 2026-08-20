//! nlr-scan integration tests: end-to-end validation of the scanning engine with the real
//! mot.txt/store.txt.

use nlr_config::MotifDefinition;
use nlr_scan::MotifParser;
use std::path::Path;

fn real_parser() -> MotifParser {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../src");
    let def = MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt")).unwrap();
    MotifParser::new(def)
}

#[test]
fn short_sequence_returns_empty() {
    let parser = real_parser();
    // Sequences shorter than the maximum motif length (50) return empty immediately.
    let list = parser.find_motifs("short", "ACDEFGHIKLMNPQRSTVWY");
    assert_eq!(list.motifs.len(), 0);
}

#[test]
fn motif_1_consensus_hits() {
    // motif_1 (P-loop) consensus = PIWGMGGVGKTTLARAVYNDP (21 aa).
    // Repeat this sequence to build length > 50; expect at least motif_1 to be hit.
    let consensus = "PIWGMGGVGKTTLARAVYNDP";
    let protein = consensus.repeat(3); // 63 aa > 50
    let parser = real_parser();
    let list = parser.find_motifs("consensus", &protein);
    // At least motif_1 should be hit.
    assert!(
        list.motifs.iter().any(|m| m.id == 1),
        "expected motif_1 to be hit, got {} motifs",
        list.motifs.len()
    );
}

#[test]
fn no_crash_on_all_stop() {
    let parser = real_parser();
    let protein = "*".repeat(60);
    let list = parser.find_motifs("stops", &protein);
    // Stop-codon sequences must not crash; may return empty or a few hits.
    let _ = list.motifs.len();
}
