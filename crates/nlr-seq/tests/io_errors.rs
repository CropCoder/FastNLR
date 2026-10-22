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
