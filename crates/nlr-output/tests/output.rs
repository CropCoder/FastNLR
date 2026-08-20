//! nlr-output integration tests: validate output formats.

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;
use nlr_core::strand::Strand;
use nlr_output::{write_nlr_bed, write_nlr_gff, write_report_txt};

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

#[test]
fn report_txt_format() {
    let nlrs = vec![sample_nlr()];
    let mut buf = Vec::new();
    write_report_txt(&mut buf, &nlrs, &def()).unwrap();
    let s = String::from_utf8(buf).unwrap();
    let cols: Vec<&str> = s.trim_end().split('\t').collect();
    assert_eq!(cols.len(), 7);
    assert_eq!(cols[0], "chr1");
    assert_eq!(cols[1], "chr1_nlr1");
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
