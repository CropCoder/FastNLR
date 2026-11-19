//! FASTA reader (streaming, gzip support).
//!
//! - gzip magic-number sniffing (first 2 bytes `0x1f 0x8b`);
//! - `read_entry` returns sequences one by one; identifier is the first whitespace-delimited token of the header with `>` stripped, description is retained;
//! - each line is `.trim()`-ed of leading/trailing whitespace before being appended to the sequence.
//! - empty/duplicate identifiers, empty records and non-sequence bytes are rejected.

use crate::translate::BioSequence;
use std::collections::HashSet;
use std::io::{self, BufRead, BufReader};

/// Parse a complete FASTA header line (`>id description`) into `(id, description)`.
///
/// The identifier is the first whitespace-delimited token after `>`; a header
/// without an identifier is rejected instead of silently becoming the empty id.
fn parse_header(line: &str) -> io::Result<(String, String)> {
    let line = line.trim_end_matches(['\r', '\n']);
    let rest = line
        .strip_prefix('>')
        .ok_or_else(|| invalid_data(format!("FASTA header does not start with '>': {line}")))?;
    parse_header_fields(rest)
}

/// Parse the description portion after `>` (shared with the streaming chopper).
pub(crate) fn parse_header_fields(rest: &str) -> io::Result<(String, String)> {
    let mut parts = rest.split_whitespace();
    let identifier = parts.next().unwrap_or("");
    if identifier.is_empty() {
        return Err(invalid_data("FASTA header has an empty identifier"));
    }
    let description = parts.collect::<Vec<_>>().join(" ");
    Ok((identifier.to_string(), description))
}

/// Whether a byte can occur in a sequence accepted by FastNLR.
///
/// Letters are accepted (DNA IUPAC and protein residues), plus `*` (stops in
/// protein input) and `-` (alignment gaps). Whitespace and non-ASCII bytes are
/// rejected so byte-based coordinates and slicing cannot get out of sync.
#[inline]
pub(crate) fn is_valid_sequence_byte(b: u8) -> bool {
    b.is_ascii_alphabetic() || b == b'*' || b == b'-'
}

pub(crate) fn invalid_data(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Streaming FASTA reader.
pub struct FastaReader<R: BufRead> {
    reader: R,
    last_line: Option<String>,
    seen_ids: HashSet<String>,
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
        Self::new(reader)
    }
}

impl<R: BufRead> FastaReader<R> {
    /// Construct from any BufRead, skipping leading non-header lines.
    pub fn new(mut reader: R) -> std::io::Result<Self> {
        let mut last_line = None;
        loop {
            let mut line = String::new();
            let bytes = reader.read_line(&mut line)?;
            if bytes == 0 {
                break; // EOF
            }
            if line.trim_start().starts_with('>') {
                last_line = Some(line);
                break;
            }
        }
        Ok(FastaReader {
            reader,
            last_line,
            seen_ids: HashSet::new(),
        })
    }

    /// Read the next sequence; returns `Ok(None)` at EOF.
    pub fn read_entry(&mut self) -> std::io::Result<Option<BioSequence>> {
        let Some(header) = self.last_line.take() else {
            return Ok(None);
        };
        let (identifier, description) = parse_header(header.trim_start())?;
        if !self.seen_ids.insert(identifier.clone()) {
            return Err(invalid_data(format!(
                "duplicate FASTA identifier '{identifier}'"
            )));
        }

        let mut sequence = String::new();
        loop {
            let mut line = String::new();
            let bytes = self.reader.read_line(&mut line)?;
            if bytes == 0 {
                break; // EOF
            }
            let trimmed = line.trim();
            if trimmed.starts_with('>') {
                self.last_line = Some(line);
                break;
            }
            validate_sequence_line(trimmed, &identifier)?;
            sequence.push_str(trimmed);
        }
        if sequence.is_empty() {
            return Err(invalid_data(format!(
                "FASTA record '{identifier}' has an empty sequence"
            )));
        }

        let mut seq = BioSequence::new(identifier, sequence);
        seq.description = description;
        Ok(Some(seq))
    }
}

/// Convenience function: read an entire FASTA file into a list of sequences.
pub fn read_all(path: &std::path::Path) -> std::io::Result<Vec<BioSequence>> {
    let mut reader = FastaReader::from_file(path)?;
    let mut out = Vec::new();
    while let Some(seq) = reader.read_entry()? {
        out.push(seq);
    }
    Ok(out)
}

fn validate_sequence_line(line: &str, identifier: &str) -> std::io::Result<()> {
    if let Some(byte) = line.bytes().find(|byte| !is_valid_sequence_byte(*byte)) {
        return Err(invalid_data(format!(
            "FASTA record '{identifier}' contains invalid sequence byte 0x{byte:02X}; \
             only ASCII letters, '*' and '-' are accepted"
        )));
    }
    Ok(())
}
