//! FASTA/chopper 输入错误传播回归测试。
//!
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use std::io::{self, BufRead, Read};

struct FailingReader;

impl Read for FailingReader {
    fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
        Err(io::Error::new(io::ErrorKind::InvalidData, "read failure"))
    }
}

impl BufRead for FailingReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        Err(io::Error::new(io::ErrorKind::InvalidData, "read failure"))
    }

    fn consume(&mut self, _amt: usize) {}
}

#[test]
fn fasta_reader_propagates_io_errors() {
    let reader: Box<dyn BufRead> = Box::new(FailingReader);
    assert!(nlr_seq::fasta::FastaReader::new(reader).is_err());
}

#[test]
fn sequence_chopper_propagates_io_errors() {
    let reader: Box<dyn BufRead> = Box::new(FailingReader);
    assert!(nlr_seq::SequenceChopper::new(reader, 10, 2).is_err());
}

#[test]
fn sequence_chopper_rejects_invalid_windows() {
    let reader: Box<dyn BufRead> = Box::new(io::Cursor::new(b">chr1\nACGT\n".to_vec()));
    assert!(nlr_seq::SequenceChopper::new(reader, 0, 0).is_err());
    let reader: Box<dyn BufRead> = Box::new(io::Cursor::new(b">chr1\nACGT\n".to_vec()));
    assert!(nlr_seq::SequenceChopper::new(reader, 10, 10).is_err());
}

#[test]
fn corrupt_gzip_is_reported() {
    let mut path = std::env::temp_dir();
    path.push(format!("fastnlr_bad_gzip_{}.fa.gz", std::process::id()));
    std::fs::write(&path, [0x1f, 0x8b, 0x00, 0x01]).unwrap();

    let fasta_err = nlr_seq::fasta::FastaReader::from_file(&path).err();
    assert!(fasta_err.is_some(), "corrupt gzip FASTA must fail");

    let chopper_err = nlr_seq::SequenceChopper::from_file(&path, 10, 2).err();
    assert!(chopper_err.is_some(), "corrupt gzip chopper must fail");
    let _ = std::fs::remove_file(path);
}

fn drain_chopper(mut chopper: nlr_seq::SequenceChopper) -> io::Result<usize> {
    let mut count = 0;
    while chopper.next_sequence()?.is_some() {
        count += 1;
    }
    Ok(count)
}

fn chopper_errors(data: &[u8]) -> (io::Error, io::Error) {
    let reader: Box<dyn BufRead> = Box::new(io::Cursor::new(data.to_vec()));
    let reader_error = match nlr_seq::SequenceChopper::new(reader, 8, 2).and_then(drain_chopper) {
        Ok(_) => panic!("reader path accepted malformed FASTA"),
        Err(error) => error,
    };

    let dir = tempfile::TempDir::new().unwrap();
    let path = dir.path().join("input.fa");
    std::fs::write(&path, data).unwrap();
    let mmap_error =
        match nlr_seq::SequenceChopper::from_mmap_file(&path, 8, 2).and_then(drain_chopper) {
            Ok(_) => panic!("mmap path accepted malformed FASTA"),
            Err(error) => error,
        };
    (reader_error, mmap_error)
}

#[test]
fn chopper_rejects_duplicate_ids_on_both_paths() {
    let (reader_error, mmap_error) = chopper_errors(b">chr1\nACGT\n>chr1\nTGCA\n");
    for error in [reader_error, mmap_error] {
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error
            .to_string()
            .contains("duplicate FASTA identifier 'chr1'"));
    }
}

#[test]
fn chopper_rejects_empty_id_on_both_paths() {
    let (reader_error, mmap_error) = chopper_errors(b">\nACGT\n");
    for error in [reader_error, mmap_error] {
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("empty identifier"));
    }
}

#[test]
fn chopper_rejects_empty_sequence_on_both_paths() {
    let (reader_error, mmap_error) = chopper_errors(b">chr1\n>chr2\nACGT\n");
    for error in [reader_error, mmap_error] {
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error
            .to_string()
            .contains("FASTA record 'chr1' has an empty sequence"));
    }
}

#[test]
fn chopper_rejects_invalid_sequence_bytes_on_both_paths() {
    let (reader_error, mmap_error) = chopper_errors(b">chr1\nAC GT\n");
    for error in [reader_error, mmap_error] {
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(error.to_string().contains("invalid sequence byte 0x20"));
    }
}

#[test]
fn chopper_rejects_non_utf8_headers_on_both_paths() {
    let (reader_error, mmap_error) = chopper_errors(b">chr\xff\nACGT\n");
    for error in [reader_error, mmap_error] {
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }
}

#[test]
fn fasta_reader_rejects_strictly_invalid_records() {
    let cases: [(&[u8], &str); 4] = [
        (b">chr1\nACGT\n>chr1\nTGCA\n", "duplicate FASTA identifier"),
        (b">\nACGT\n", "empty identifier"),
        (b">chr1\n>chr2\nACGT\n", "has an empty sequence"),
        (b">chr1\nAC GT\n", "invalid sequence byte 0x20"),
    ];
    for (data, expected) in cases {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("input.fa");
        std::fs::write(&path, data).unwrap();
        let error = nlr_seq::fasta::read_all(&path)
            .expect_err("strict FASTA reader must reject malformed input");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(
            error.to_string().contains(expected),
            "expected {expected:?}, got {error}"
        );
    }
}
