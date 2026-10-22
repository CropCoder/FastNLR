//! nlr-assemble integration tests: validate the three-step assembly.

use nlr_assemble::{assemble, AssembleParams};
use nlr_core::motif::Motif;
use nlr_core::signature_def::{AnnotatorSignatureDefinition, FlexibleSeedConfig};
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
fn rnl_cc_ploop_seed_recovers_complete_locus() {
    // ADR1-L1 形态：helper CC motif 位于内部命中串中，有效的精确 seed 是 [23,1]。
    // 后续 LRR 通过 elongation 并入，最终应输出一个完整 CC-NBARC-LRR。
    let motifs = vec![
        dna_motif(21, 1000, Strand::Forward),
        dna_motif(23, 1100, Strand::Forward),
        dna_motif(1, 1500, Strand::Forward),
        dna_motif(11, 2000, Strand::Forward),
        dna_motif(11, 2100, Strand::Forward),
        dna_motif(11, 2200, Strand::Forward),
        dna_motif(11, 2300, Strand::Forward),
    ];
    let d = def();
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &d);
    assert_eq!(nlrs.len(), 1);
    assert_eq!(nlrs[0].domain_string(&d), "CC-NBARC-LRR");
}

#[test]
fn rnl_cc_without_ploop_does_not_seed() {
    // 精度保护：新增补丁不能仅凭 helper CC 或“CC+LRR”造位点。
    let motifs = vec![
        dna_motif(21, 1000, Strand::Forward),
        dna_motif(23, 1100, Strand::Forward),
        dna_motif(11, 1300, Strand::Forward),
    ];
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &def());
    assert!(nlrs.is_empty());
}

#[test]
fn unresolved_cc_ploop_fragment_is_suppressed() {
    // 单个 CC 与 P-loop 可以形成 seed，但若没有 LRR/第三 motif，实际只是
    // 200-500 bp 的片段；默认输出不应将其报告为一个独立 NLR 位点。
    let motifs = vec![
        dna_motif(21, 1000, Strand::Forward),
        dna_motif(1, 1200, Strand::Forward),
    ];
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &def());
    assert!(nlrs.is_empty());
}

#[test]
fn cc_ploop_with_third_motif_survives_as_partial_locus() {
    // DAR5 形态：两个 CC motif + P-loop，虽无 LRR，但已超出无信息的两 motif 片段。
    let motifs = vec![
        dna_motif(21, 1000, Strand::Forward),
        dna_motif(23, 1100, Strand::Forward),
        dna_motif(1, 1500, Strand::Forward),
    ];
    let d = def();
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &d);
    assert_eq!(nlrs.len(), 1);
    assert_eq!(nlrs[0].domain_string(&d), "CC-NBARC");
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

#[test]
fn flexible_seed_is_additive_not_replacing() {
    // 回归测试（2026-10-12）：开启 --flexible-seed 后必须**同时**保留精确种子组合的播种能力。
    //
    // 构造 RefPlantNLR 中 Pb1 的形状：命中精确组合 [12,2,8]，但其中 motif_8 属 linker、
    // 真正的 NB-ARC 类 motif 只有 motif_12 与 motif_2 共 2 个，达不到柔性门槛 min_nbarc=3。
    // 修复前播种是 if/else 链——开了柔性就整体顶替精确组合，这条命中串反而播不了种，
    // 出现「默认能检出、开 --flexible-seed 反而漏检」的回退。
    let motifs = vec![
        dna_motif(12, 1000, Strand::Forward),
        dna_motif(2, 1100, Strand::Forward),
        dna_motif(8, 1200, Strand::Forward),
    ];
    let params = AssembleParams::default();

    // 默认（精确组合）能播种，作为基准
    let plain = assemble("chr1", motifs.clone(), &params, &def());
    assert_eq!(plain.len(), 1, "精确组合 [12,2,8] 应能播种");

    // 开启柔性播种后，同一命中串仍应播种（并集语义）
    let flex_def = def().with_flexible_seed(FlexibleSeedConfig::default());
    let flex = assemble("chr1", motifs, &params, &flex_def);
    assert_eq!(
        flex.len(),
        1,
        "开启 --flexible-seed 不应丢失精确组合能检出的位点"
    );
    let ids: Vec<u8> = flex[0].motifs.iter().map(|m| m.id).collect();
    assert_eq!(ids, vec![12, 2, 8]);
}

#[test]
fn flexible_seed_still_recovers_rank_order_loci() {
    // 反向保护：并集不得削弱柔性播种本身的能力。
    // [6,4,10,3,2] 的 rank 为 5,6,8,9,11 单调递增、NB-ARC 数 5 >= 3，但不匹配任何精确组合，
    // 因此只有柔性条件能检出。
    let motifs = vec![
        dna_motif(6, 1000, Strand::Forward),
        dna_motif(4, 1100, Strand::Forward),
        dna_motif(10, 1200, Strand::Forward),
        dna_motif(3, 1300, Strand::Forward),
        dna_motif(2, 1400, Strand::Forward),
    ];
    let params = AssembleParams::default();

    let plain = assemble("chr1", motifs.clone(), &params, &def());
    assert!(plain.is_empty(), "该串不匹配精确组合，默认配置不应播种");

    let flex_def = def().with_flexible_seed(FlexibleSeedConfig::default());
    let flex = assemble("chr1", motifs, &params, &flex_def);
    assert_eq!(flex.len(), 1, "柔性播种应能检出 rank 单调的 NB-ARC 串");
}

#[test]
fn flexible_seed_one_accepts_an_isolated_ploop() {
    let motifs = vec![dna_motif(1, 1000, Strand::Forward)];
    let d = def().with_flexible_seed(FlexibleSeedConfig {
        min_nbarc: 1,
        require_ploop: true,
    });
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &d);
    assert_eq!(nlrs.len(), 1, "single P-loop should seed with min_nbarc=1");
    assert_eq!(nlrs[0].motifs.len(), 1);
}

#[test]
fn absorbs_nested_tandem_repeat_unit() {
    // RPP2A 形态（金标准真值 TIR-NBARC-TIR-NBARC-LRR）：串联重复的第二单元因 rank 回退
    // 无法被 merge_seeds 合并，会自成一个小位点并**嵌套在主位点跨度内**。
    // 组装后处理应把嵌套位点吸收回主位点，最终只输出一个位点且域串完整。
    let motifs = vec![
        // 主位点：TIR + NB-ARC（前半段）……
        dna_motif(18, 1000, Strand::Forward),
        dna_motif(15, 1100, Strand::Forward),
        dna_motif(13, 1200, Strand::Forward),
        dna_motif(1, 1300, Strand::Forward),
        dna_motif(4, 1400, Strand::Forward),
        dna_motif(5, 1500, Strand::Forward),
        dna_motif(10, 1600, Strand::Forward),
        dna_motif(3, 1700, Strand::Forward),
        // ……以及 LRR（后半段），中间留出的区间正好容纳被嵌套的第二单元。
        dna_motif(11, 2900, Strand::Forward),
        dna_motif(11, 3000, Strand::Forward),
        dna_motif(9, 3100, Strand::Forward),
        dna_motif(11, 3200, Strand::Forward),
        // 第二单元：TIR + NB-ARC，完全落在主位点跨度 (1000, 3209) 之内。
        dna_motif(15, 2000, Strand::Forward),
        dna_motif(13, 2100, Strand::Forward),
        dna_motif(1, 2200, Strand::Forward),
        dna_motif(4, 2300, Strand::Forward),
        dna_motif(5, 2400, Strand::Forward),
        dna_motif(10, 2500, Strand::Forward),
        dna_motif(3, 2600, Strand::Forward),
    ];
    let d = def();
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &d);
    assert_eq!(nlrs.len(), 1, "嵌套的重复单元应被吸收，只保留一个位点");
    assert_eq!(nlrs[0].domain_string(&d), "TIR-NBARC-TIR-NBARC-LRR");
}

#[test]
fn absorbs_trailing_nbarc_segment() {
    // Prv 形态（金标准真值 TIR-NBARC-LRR-NBARC）：尾随 NB-ARC 的 rank(4) 低于 LRR(14)，
    // elongate 与 merge_seeds 都不会把它并入，故形成下游孤立的 NB-ARC 段。
    // 组装后处理应把它按串联结构补回。
    let motifs = vec![
        dna_motif(18, 1000, Strand::Forward),
        dna_motif(15, 1100, Strand::Forward),
        dna_motif(1, 1200, Strand::Forward),
        dna_motif(4, 1300, Strand::Forward),
        dna_motif(5, 1400, Strand::Forward),
        dna_motif(10, 1500, Strand::Forward),
        dna_motif(3, 1600, Strand::Forward),
        dna_motif(12, 1700, Strand::Forward),
        dna_motif(7, 1800, Strand::Forward),
        dna_motif(11, 1900, Strand::Forward),
        dna_motif(11, 2000, Strand::Forward),
        dna_motif(9, 2100, Strand::Forward),
        dna_motif(11, 2200, Strand::Forward),
        // 下游 180 bp 处的孤立 NB-ARC 段（<= distance_for_elongating = 2500）
        dna_motif(1, 2400, Strand::Forward),
        dna_motif(4, 2500, Strand::Forward),
        dna_motif(5, 2600, Strand::Forward),
        dna_motif(10, 2700, Strand::Forward),
        dna_motif(3, 2800, Strand::Forward),
    ];
    let d = def();
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &d);
    assert_eq!(nlrs.len(), 1, "尾随 NB-ARC 段应被并入主位点");
    assert_eq!(nlrs[0].domain_string(&d), "TIR-NBARC-LRR-NBARC");
}

#[test]
fn absorbs_trailing_nbarc_segment_on_reverse_strand() {
    let motifs = vec![
        dna_motif(18, 3000, Strand::Reverse),
        dna_motif(15, 2900, Strand::Reverse),
        dna_motif(1, 2700, Strand::Reverse),
        dna_motif(4, 2600, Strand::Reverse),
        dna_motif(5, 2500, Strand::Reverse),
        dna_motif(10, 2400, Strand::Reverse),
        dna_motif(3, 2300, Strand::Reverse),
        dna_motif(12, 2200, Strand::Reverse),
        dna_motif(7, 2100, Strand::Reverse),
        dna_motif(11, 2000, Strand::Reverse),
        dna_motif(11, 1900, Strand::Reverse),
        dna_motif(9, 1800, Strand::Reverse),
        dna_motif(11, 1700, Strand::Reverse),
        // Reverse-strand C-terminal tail is on the lower-coordinate side.
        dna_motif(1, 1500, Strand::Reverse),
        dna_motif(4, 1400, Strand::Reverse),
        dna_motif(5, 1300, Strand::Reverse),
        dna_motif(10, 1200, Strand::Reverse),
        dna_motif(3, 1100, Strand::Reverse),
    ];
    let d = def();
    let nlrs = assemble("chr1", motifs, &AssembleParams::default(), &d);
    assert_eq!(nlrs.len(), 1, "reverse C-terminal tail should merge");
    assert_eq!(nlrs[0].domain_string(&d), "TIR-NBARC-LRR-NBARC");
}
