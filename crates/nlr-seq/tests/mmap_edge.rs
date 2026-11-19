//! nlr-seq mmap chopper edge-case test: verifies the mmap path and reader path produce consistent output for various inputs.

use nlr_seq::SequenceChopper;
use std::io::Write;

/// (序列名, 序列) 列表，mmap 与 reader 两条路径的结果类型。
type Fragments = Vec<(String, String)>;

/// Chop the same input via both mmap and reader paths, returning both results.
fn both_paths(data: &[u8], frag: usize, overlap: usize) -> (Fragments, Fragments) {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let uniq = COUNTER.fetch_add(1, Ordering::SeqCst);

    // mmap path (unique filename to avoid concurrency conflicts).
    let mut path = std::env::temp_dir();
    path.push(format!("nlr_edge_{}_{}.fasta", std::process::id(), uniq));
    {
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(data).unwrap();
    }
    let mut c_mmap = SequenceChopper::from_file(&path, frag, overlap).unwrap();
    let mut mmap_out = Vec::new();
    while let Some(s) = c_mmap.next_sequence().unwrap() {
        mmap_out.push((s.identifier.clone(), s.sequence.clone()));
    }
    let _ = std::fs::remove_file(&path);

    // reader path.
    let reader: Box<dyn std::io::BufRead> = Box::new(std::io::Cursor::new(data.to_vec()));
    let mut c_reader = SequenceChopper::new(reader, frag, overlap).unwrap();
    let mut reader_out = Vec::new();
    while let Some(s) = c_reader.next_sequence().unwrap() {
        reader_out.push((s.identifier.clone(), s.sequence.clone()));
    }

    (mmap_out, reader_out)
}

fn assert_consistent(name: &str, data: &[u8], frag: usize, overlap: usize) {
    let (mmap_out, reader_out) = both_paths(data, frag, overlap);
    assert_eq!(
        mmap_out, reader_out,
        "[{}] mmap and reader paths inconsistent\nmmap:   {:?}\nreader: {:?}",
        name, mmap_out, reader_out
    );
}

fn assert_rejected_consistently(data: &[u8], expected: &str) {
    use std::sync::atomic::{AtomicU64, Ordering};
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    let uniq = COUNTER.fetch_add(1, Ordering::SeqCst);

    let mut path = std::env::temp_dir();
    path.push(format!("nlr_reject_{}_{}.fasta", std::process::id(), uniq));
    {
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(data).unwrap();
    }
    let mmap_error = match SequenceChopper::from_file(&path, 10, 4)
        .and_then(|mut chopper| chopper.next_sequence().map(|entry| entry.map(|_| ())))
    {
        Ok(_) => panic!("mmap path accepted invalid FASTA"),
        Err(error) => error,
    };
    let _ = std::fs::remove_file(&path);

    let reader: Box<dyn std::io::BufRead> = Box::new(std::io::Cursor::new(data.to_vec()));
    let reader_error = match SequenceChopper::new(reader, 10, 4)
        .and_then(|mut chopper| chopper.next_sequence().map(|entry| entry.map(|_| ())))
    {
        Ok(_) => panic!("reader path accepted invalid FASTA"),
        Err(error) => error,
    };

    for error in [mmap_error, reader_error] {
        assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
        assert!(error.to_string().contains(expected));
    }
}

#[test]
fn empty_file_consistent() {
    assert_consistent("empty", b"", 10, 4);
}

#[test]
fn header_only_is_rejected_consistently() {
    assert_rejected_consistently(b">chr1", "empty sequence");
    assert_rejected_consistently(b">chr1\n", "empty sequence");
}

#[test]
fn windows_crlf_consistent() {
    // CRLF line endings: \r should be skipped, not mixed into the sequence.
    assert_consistent("crlf", b">chr1\r\nACGTACGTACGT\r\n", 4, 2);
}

#[test]
fn no_trailing_newline_consistent() {
    assert_consistent("notail", b">chr1\nACGTACGT", 4, 2);
}

#[test]
fn multi_sequence_consistent() {
    assert_consistent("multi", b">chr1\nACGTACGTACGT\n>chr2\nTTTTGGGG\n", 4, 2);
}

#[test]
fn crlf_does_not_pollute_sequence() {
    // Under CRLF the first fragment should be exactly ACGT (no \r).
    let (mmap_out, _) = both_paths(b">chr1\r\nACGTACGTACGT\r\n", 4, 2);
    assert_eq!(
        mmap_out[0].1, "ACGT",
        "first fragment should not contain \\r, got {:?}",
        mmap_out[0].1
    );
    // 12 bases, frag=4 overlap=2 -> 5 fragments.
    assert_eq!(mmap_out.len(), 5);
    // No fragment should contain \r.
    for (_, seq) in &mmap_out {
        assert!(!seq.contains('\r'), "sequence contains \\r: {:?}", seq);
    }
}
