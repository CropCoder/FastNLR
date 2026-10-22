//! nlr-output integration tests: validate output formats.

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;
use nlr_core::strand::Strand;
use nlr_output::{
    export_motifs, import_motifs, write_nbarc_alignment_fasta, write_nlr_bed, write_nlr_gff,
    write_report_txt,
};

fn def() -> AnnotatorSignatureDefinition {
    AnnotatorSignatureDefinition::new()
}

fn sample_nlr() -> MotifList {
    // A complete NLR: P-loop(1) + NBARC(6) + LRR(9), forward strand.
    let mut m1 = Motif::new_protein(1, "chr1".to_string(), 1, "P".to_string(), 1e-6);
    m1.set_dna("chr1".to_string(), 100, 20000, 0, Strand::Forward);
    let mut m2 = Motif::new_protein(6, "chr1".to_string(), 10, "N".to_string(), 1e-6);
    m2.set_dna("chr1".to_string(), 200, 20000, 0, Strand::Forward);
    let mut m3 = Motif::new_protein(9, "chr1".to_string(), 20, "L".to_string(), 1e-6);
    m3.set_dna("chr1".to_string(), 300, 20000, 0, Strand::Forward);
    MotifList::new("chr1_nlr1".to_string(), vec![m1, m2, m3])
}

/// 构造一个宽度等于该 motif 共识长度的命中，便于按共识长度核算比对长度。
/// LRR 等无共识序列的 motif 退化为宽度 1（不参与 NB-ARC 比对长度）。
fn nbarc_motif(def: &AnnotatorSignatureDefinition, id: u8, start: u64) -> Motif {
    let width = def.consensus(id).len().max(1);
    let mut m = Motif::new_protein(id, "chr1".to_string(), 1, "A".repeat(width), 1e-6);
    m.set_dna("chr1".to_string(), start, 20000, 0, Strand::Forward);
    m
}

/// 8 个 NB-ARC motif 共识长度之和（21+29+15+20+15+21+15+29）。
const NBARC_ALIGNMENT_LEN: usize = 165;

#[test]
fn report_txt_format() {
    let nlrs = vec![sample_nlr()];
    let mut buf = Vec::new();
    write_report_txt(&mut buf, &nlrs, &def()).unwrap();
    let s = String::from_utf8(buf).unwrap();
    let cols: Vec<&str> = s.trim_end().split('\t').collect();
    // 自 2026-10-12 起不再输出 motif 编号清单：固定 6 列。
    assert_eq!(cols.len(), 6, "nlr.txt 应不含 motif 编号清单: {cols:?}");
    assert_eq!(cols[0], "chr1");
    assert_eq!(cols[1], "chr1_nlr1");
    assert_eq!(cols[5], "+");
}

#[test]
fn gff_start_is_1_based() {
    let nlrs = vec![sample_nlr()];
    let mut buf = Vec::new();
    write_nlr_gff(&mut buf, &nlrs, &def(), "2026-01-01 00:00", false).unwrap();
    let s = String::from_utf8(buf).unwrap();
    assert!(s.starts_with("##gff-version 2"));
    let body = s.lines().find(|l| l.starts_with("chr1\t")).unwrap();
    let cols: Vec<&str> = body.split('\t').collect();
    // The start column should be 1-based (101), feature column NBSLRR.
    assert_eq!(cols[2], "NBSLRR");
    assert_eq!(cols[3], "101"); // span start=100 -> +1
}

#[test]
fn bed_has_12_columns() {
    let nlrs = vec![sample_nlr()];
    let mut buf = Vec::new();
    write_nlr_bed(&mut buf, &nlrs, &def()).unwrap();
    let s = String::from_utf8(buf).unwrap();
    let body = s.lines().find(|l| l.starts_with("chr1\t")).unwrap();
    let cols: Vec<&str> = body.split('\t').collect();
    assert_eq!(cols.len(), 12);
    // A complete NLR should be green.
    assert_eq!(cols[8], "0,255,0");
}

#[test]
fn nbarc_alignment_survives_duplicate_motifs() {
    let d = def();
    // 含一条完全重复的 motif_6（同坐标），以及夹在 NB-ARC 块之间的 LRR(motif_9)。
    let spec: [(u8, u64); 10] = [
        (1, 100),
        (6, 200),
        (6, 200),
        (4, 300),
        (5, 400),
        (10, 500),
        (9, 600),
        (3, 700),
        (12, 800),
        (2, 900),
    ];
    let motifs: Vec<Motif> = spec.iter().map(|&(id, s)| nbarc_motif(&d, id, s)).collect();
    let nlr = MotifList::new("chr1_nlr1".to_string(), motifs);

    let mut buf = Vec::new();
    write_nbarc_alignment_fasta(&mut buf, &[nlr], &d, false).unwrap();
    let s = String::from_utf8(buf).unwrap();
    let seq = s.lines().nth(1).unwrap().trim_end();

    assert_eq!(seq.len(), NBARC_ALIGNMENT_LEN);
    assert_eq!(
        seq.matches('-').count(),
        0,
        "重复 motif 不应导致后续块全部变成 gap"
    );
}

#[test]
fn nbarc_alignment_pads_missing_block() {
    let d = def();
    // 缺 motif_4（共识长度 15），其余 7 个 NB-ARC 块齐全。
    let spec: [(u8, u64); 8] = [
        (1, 100),
        (6, 200),
        (5, 300),
        (10, 400),
        (3, 500),
        (12, 600),
        (2, 700),
        (9, 800),
    ];
    let motifs: Vec<Motif> = spec.iter().map(|&(id, s)| nbarc_motif(&d, id, s)).collect();
    let nlr = MotifList::new("chr1_nlr1".to_string(), motifs);

    let mut buf = Vec::new();
    write_nbarc_alignment_fasta(&mut buf, &[nlr], &d, false).unwrap();
    let s = String::from_utf8(buf).unwrap();
    let seq = s.lines().nth(1).unwrap().trim_end();

    assert_eq!(seq.len(), NBARC_ALIGNMENT_LEN);
    assert_eq!(seq.matches('-').count(), d.consensus(4).len());
}

#[test]
fn motif_checkpoint_import_is_strict() {
    let motifs = sample_nlr().motifs;
    let mut buf = Vec::new();
    export_motifs(&mut buf, &motifs).unwrap();
    let text = String::from_utf8(buf).unwrap();
    let imported = import_motifs(&text).unwrap();
    assert_eq!(imported.len(), motifs.len());
    assert!(imported.iter().all(|m| m.dna_parameters_set));

    assert!(import_motifs("not-a-checkpoint-row\n").is_err());

    let no_dna = Motif::new_protein(1, "protein".to_string(), 1, "P".to_string(), 1e-6);
    let mut buf = Vec::new();
    export_motifs(&mut buf, &[no_dna]).unwrap();
    let text = String::from_utf8(buf).unwrap();
    assert!(
        import_motifs(&text).is_err(),
        "DNA-free motif must not enter a genome checkpoint"
    );
}
