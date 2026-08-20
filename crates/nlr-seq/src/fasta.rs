//! FASTA reader (streaming, gzip support).
//!
//! Reimplementation of Java `FastaReader`:
//! - gzip magic-number sniffing (first 2 bytes `0x1f 0x8b`);
//! - `read_entry` returns sequences one by one; identifier is the first whitespace-delimited token of the header with `>` stripped, description is retained;
//! - each line is `.trim()`-ed of leading/trailing whitespace before being appended to the sequence.

use crate::translate::BioSequence;
use std::io::{BufRead, BufReader};

/// Streaming FASTA reader.
pub struct FastaReader<R: BufRead> {
    reader: R,
    last_line: Option<String>,
}

impl FastaReader<Box<dyn BufRead>> {
    /// Construct from file, auto-sniffing gzip.
    pub fn from_file(path: &std::path::Path) -> std::io::Result<Self> {
        let mut file = std::fs::File::open(path)?;
        let mut header = [0u8; 2];
        let is_gzip = {
            use std::io::Read;
            let n = file.read(&mut header)?;
            n == 2 && header == [0x1f, 0x8b]
        };
        // Reset offset after reading magic bytes (try_clone shares offset on Unix; use seek back to 0 instead).
        use std::io::Seek;
        file.seek(std::io::SeekFrom::Start(0))?;
        let reader: Box<dyn BufRead> = if is_gzip {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };
        Ok(Self::new(reader))
    }
}

impl<R: BufRead> FastaReader<R> {
    /// Construct from any BufRead, skipping leading non-header lines.
    pub fn new(mut reader: R) -> Self {
        let mut last_line = None;
        loop {
            let mut line = String::new();
            match reader.read_line(&mut line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    if line.trim_start().starts_with('>') {
                        last_line = Some(line);
                        break;
                    }
                }
                Err(_) => break,
            }
        }
        FastaReader { reader, last_line }
    }

    /// Read the next sequence; returns None at EOF.
    pub fn read_entry(&mut self) -> Option<BioSequence> {
        let header = self.last_line.take()?;
        let header = header.trim_end();
        // identifier = first whitespace-delimited token with '>' stripped; description = remainder.
        let mut parts = header.split_whitespace();
        let identifier = parts.next()?.trim_start_matches('>').to_string();
        let description = parts.collect::<Vec<_>>().join(" ");

        let mut sequence = String::new();
        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    let trimmed = line.trim();
                    if trimmed.starts_with('>') {
                        self.last_line = Some(line);
                        break;
                    }
                    sequence.push_str(trimmed);
                }
                Err(_) => break,
            }
        }

        let mut seq = BioSequence::new(identifier, sequence);
        seq.description = description;
        Some(seq)
    }
}

/// Convenience function: read an entire FASTA file into a list of sequences.
pub fn read_all(path: &std::path::Path) -> std::io::Result<Vec<BioSequence>> {
    let mut reader = FastaReader::from_file(path)?;
    let mut out = Vec::new();
    while let Some(seq) = reader.read_entry() {
        out.push(seq);
    }
    Ok(out)
}
