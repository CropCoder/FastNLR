//! nlr-seq mmap chopper integration test: verifies the mmap path and reader path produce identical output.

use nlr_seq::SequenceChopper;
use std::io::Write;

#[test]
fn mmap_chopper_matches_reader() {
    // Build a multi-line, multi-sequence FASTA.
    let mut fasta = Vec::new();
    fasta.extend_from_slice(b">chr1\n");
    for _ in 0..3 {
        fasta.extend_from_slice(b"ACGTACGTACGTACGTACGTACGT\n");
    }
    fasta.extend_from_slice(b">chr2\n");
    fasta.extend_from_slice(b"TTTTGGGGCCCCAAAATTTTGGGG\n");

    // Write to a temporary file.
    let mut path = std::env::temp_dir();
    path.push(format!("nlr_mmap_test_{}.fasta", std::process::id()));
    {
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(&fasta).unwrap();
    }

    // mmap path (from_file auto-selects mmap for non-gzip).
    let mut c = SequenceChopper::from_file(&path, 10, 4).unwrap();
    let mut mmap_seqs = Vec::new();
    while let Some(s) = c.next_sequence() {
        mmap_seqs.push((s.identifier.clone(), s.sequence.clone()));
    }

    // reader path (constructed directly with a Cursor).
    let reader: Box<dyn std::io::BufRead> = Box::new(std::io::Cursor::new(fasta.clone()));
    let mut c2 = SequenceChopper::new(reader, 10, 4).unwrap();
    let mut reader_seqs = Vec::new();
    while let Some(s) = c2.next_sequence() {
        reader_seqs.push((s.identifier.clone(), s.sequence.clone()));
    }

    // Both paths should produce identical output.
    assert_eq!(mmap_seqs, reader_seqs, "mmap and reader paths produced inconsistent output");

    let _ = std::fs::remove_file(&path);
}
