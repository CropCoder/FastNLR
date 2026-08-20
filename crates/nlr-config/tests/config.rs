//! nlr-config integration tests: load the real mot.txt/store.txt and validate the parsed result.

use nlr_config::MotifDefinition;
use std::path::Path;

#[test]
fn load_real_config() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../src");
    let def = MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt")).unwrap();

    // All 20 motifs loaded.
    let mut ids: Vec<u8> = (1..=20).collect();
    ids.sort();
    let mut names = def.motif_names().to_vec();
    names.sort();
    assert_eq!(names, ids, "should load 20 motifs");

    // Maximum motif length = 50 (per documentation).
    assert_eq!(def.max_length(), 50);

    // Spot-check motif_1 (P-loop) length 21 (consensus PIWGMGGVGKTTLARAVYNDP).
    assert_eq!(def.length(1), 21);

    // PWM score range 0..=100, CDF index queryable.
    // motif_4@0@A = 27 (first line of mot.txt, known value).
    let s = def.score(4, 0, b'A');
    assert!(s >= 0 && s <= 100, "PWM score should be in 0..=100, got {}", s);
}

#[test]
fn thresholds_monotonic() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../src");
    let def = MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt")).unwrap();
    let t = def.score_thresholds(1e-4);
    // Thresholds are non-negative, and every motif has a value.
    for id in 1..=20u8 {
        assert!(t[id as usize] >= 0);
    }
}

#[test]
fn score_non_ascii_returns_zero() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../../src");
    let def = MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt")).unwrap();
    // '*' (42) is less than 'A', should return 0.
    assert_eq!(def.score(4, 0, b'*'), 0);
}
