//! nlr-config integration tests: load the real mot.txt/store.txt and validate the parsed result.
//!
//! NOTE: the data files live in `crates/nlr-cli/data/` (the crate layout moved there long ago);
//! the previous path (`../../../src`) pointed outside the repository, so these three tests failed
//! on every CI run. The path below is relative to this crate's manifest dir.

use nlr_config::MotifDefinition;
use std::path::Path;

#[test]
fn load_real_config() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../nlr-cli/data");
    let def = MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt")).unwrap();

    // The built-in library has the standard 20 motifs; an extended library may have more,
    // so assert "the 20 built-ins are present" instead of an exact count.
    let mut names = def.motif_names().to_vec();
    names.sort();
    let builtins: Vec<u8> = (1..=20).collect();
    assert!(names.len() >= builtins.len(), "should load at least 20 motifs, got {}", names.len());
    assert!(builtins.iter().all(|i| names.contains(i)), "missing built-in motifs");
    assert_eq!(def.max_motif_id(), *names.last().unwrap());

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
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../nlr-cli/data");
    let def = MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt")).unwrap();
    let t = def.score_thresholds(1e-4);
    // Thresholds are non-negative, and every motif has a value.
    for id in 1..=20u8 {
        assert!(t[id as usize] >= 0);
    }
}

#[test]
fn score_non_ascii_returns_zero() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../nlr-cli/data");
    let def = MotifDefinition::load(&dir.join("mot.txt"), &dir.join("store.txt")).unwrap();
    // '*' (42) is less than 'A', should return 0.
    assert_eq!(def.score(4, 0, b'*'), 0);
}

#[test]
fn missing_cdf_table_is_an_error_not_a_panic() {
    // mot.txt 声明两个 motif，store.txt 只覆盖 motif_1：旧实现会在 score_thresholds() 越界 panic。
    let mot = "motif_1@0@G 10\nmotif_2@0@K 10\n";
    let store = "motif_1@0 1.0\nmotif_1@10 1e-9\n";
    let err = match MotifDefinition::load_from_str(mot, store) {
        Ok(_) => panic!("mot/store 不一致时应返回错误"),
        Err(e) => e,
    };
    let msg = err.to_string();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    assert!(msg.contains("motif_2"), "错误信息应指明缺失的 motif: {msg}");
    assert!(msg.contains("store.txt"), "错误信息应指明问题文件: {msg}");
}

#[test]
fn cdf_out_of_range_uses_last_tabulated_value() {
    // 表只覆盖到 score=5，最后一个表项是 1e-3。
    let mot = "motif_1@0@G 10\n";
    let store = "motif_1@0 1.0\nmotif_1@5 1e-3\n";
    let def = MotifDefinition::load_from_str(mot, store).unwrap();
    // 表内查询照旧。
    assert_eq!(def.cdf(1, 5), 1e-3);
    // 越界查询返回表内最后一个（最小）p 值，而不是旧实现的 0.0。
    assert_eq!(def.cdf(1, 9999), 1e-3);
    assert_ne!(def.cdf(1, 9999), 0.0);
}

#[test]
fn score_thresholds_never_matching_marks_motif_unreachable() {
    // 表内所有 p 值都 >= 阈值 1e-4，说明该 motif 不可能达到该显著性。
    let mot = "motif_1@0@G 10\n";
    let store = "motif_1@0 1.0\nmotif_1@5 1.0\n";
    let def = MotifDefinition::load_from_str(mot, store).unwrap();
    let t = def.score_thresholds(1e-4);
    assert_eq!(t[1], i32::MAX, "无任何分数达标时应标记为不可达");
}
