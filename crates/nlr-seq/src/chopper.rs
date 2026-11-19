//! Sequence chopper (streaming slicing with overlap).
//!
//! - reads character by character, accumulating up to `fragmentLength` then chopping;
//! - fragment naming `{id}_{offset}`, where offset is the 0-based start of the fragment, step `fragmentLength - overlap`;
//! - sequence is uppercased; newlines are skipped; header takes only the first whitespace-delimited token with `>` stripped.
//! - identifiers are unique and non-empty; empty records and non-sequence bytes are rejected.
//!
//! Both the mmap and streaming paths apply the same FASTA validation.

use crate::fasta::{invalid_data, is_valid_sequence_byte, parse_header_fields};
use crate::translate::BioSequence;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};

/// Sequence chopper.
pub struct SequenceChopper {
    reader: Box<dyn BufRead>,
    fragment_length: usize,
    overlap: usize,
    current_id: String,
    current_sequence: String,
    offset: usize,
    eof: bool,
    /// mmap path: (mapping, current read position). When non-None, zero-copy slicing is used.
    mmap: Option<(memmap2::Mmap, usize)>,
    /// Identifiers already seen, to reject duplicate FASTA records.
    seen_ids: HashSet<String>,
}

impl SequenceChopper {
    fn validate_window(fragment_length: usize, overlap: usize) -> std::io::Result<()> {
        if fragment_length == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "fragment_length must be greater than 0",
            ));
        }
        if overlap >= fragment_length {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!(
                    "overlap ({overlap}) must be smaller than fragment_length ({fragment_length})"
                ),
            ));
        }
        Ok(())
    }

    /// Construct from file, auto-sniffing gzip.
    /// Non-gzip files take the mmap path (zero-copy slicing); gzip takes the streaming decompression path.
    pub fn from_file(
        path: &std::path::Path,
        fragment_length: usize,
        overlap: usize,
    ) -> std::io::Result<Self> {
        Self::validate_window(fragment_length, overlap)?;
        let mut file = std::fs::File::open(path)?;
        let mut header = [0u8; 2];
        let is_gzip = {
            use std::io::Read;
            let n = file.read(&mut header)?;
            n == 2 && header == [0x1f, 0x8b]
        };
        // Reset offset after reading magic bytes (try_clone shares offset on Unix).
        use std::io::Seek;
        file.seek(std::io::SeekFrom::Start(0))?;
        if is_gzip {
            let reader: Box<dyn BufRead> =
                Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)));
            Self::new(reader, fragment_length, overlap)
        } else {
            Self::from_mmap_file(path, fragment_length, overlap)
        }
    }

    /// mmap large-file chopper (non-gzip).
    ///
    /// Maps the whole file at once; byte-by-byte traversal incurs zero syscalls and zero copies (slices borrow; copies only on uppercase conversion as needed).
    pub fn from_mmap_file(
        path: &std::path::Path,
        fragment_length: usize,
        overlap: usize,
    ) -> std::io::Result<Self> {
        Self::validate_window(fragment_length, overlap)?;
        let file = std::fs::File::open(path)?;
        // SAFETY: read-only mapping, lifetime owned by Mmap.
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        let data = mmap.as_ref();
        let mut pos = 0usize;

        // Skip leading non-'>' content.
        while pos < data.len() && data[pos] != b'>' {
            pos += 1;
        }
        // Read identifier (from '>' to newline).
        let mut current_id = String::new();
        if pos < data.len() {
            pos += 1; // skip '>'
            let (identifier, next_pos) = parse_mmap_header(data, pos)?;
            current_id = identifier;
            pos = next_pos;
        }
        let mut seen_ids = HashSet::new();
        if !current_id.is_empty() {
            seen_ids.insert(current_id.clone());
        }

        Ok(SequenceChopper {
            reader: Box::new(std::io::empty()), // mmap path does not use reader
            fragment_length,
            overlap,
            current_id,
            current_sequence: String::new(),
            offset: 0,
            eof: false,
            mmap: Some((mmap, pos)),
            seen_ids,
        })
    }

    /// Construct from any BufRead.
    pub fn new(
        mut reader: Box<dyn BufRead>,
        fragment_length: usize,
        overlap: usize,
    ) -> std::io::Result<Self> {
        Self::validate_window(fragment_length, overlap)?;
        // Skip leading non-'>' content.
        let mut c = reader.read_byte()?;
        while let Some(ch) = c {
            if ch == b'>' {
                break;
            }
            c = reader.read_byte()?;
        }
        // Now c == Some(b'>'); read identifier.
        let current_id = if c == Some(b'>') {
            reader.read_header_id()?
        } else {
            String::new() // empty file
        };
        let mut seen_ids = HashSet::new();
        if !current_id.is_empty() {
            seen_ids.insert(current_id.clone());
        }

        Ok(SequenceChopper {
            reader,
            fragment_length,
            overlap,
            current_id,
            current_sequence: String::new(),
            offset: 0,
            eof: false,
            mmap: None,
            seen_ids,
        })
    }

    /// Read the next fragment; returns `Ok(None)` at EOF.
    pub fn next_sequence(&mut self) -> std::io::Result<Option<BioSequence>> {
        if self.eof {
            return Ok(None);
        }
        if self.mmap.is_some() {
            return self.next_sequence_mmap();
        }
        self.next_sequence_reader()
    }

    /// mmap path: bytes come from the mmap slice (saves syscalls); accumulation logic matches the reader path.
    fn next_sequence_mmap(&mut self) -> std::io::Result<Option<BioSequence>> {
        // Take out the current mapping and position.
        let Some(mmap) = self.mmap.take() else {
            return Ok(None);
        };
        let data = mmap.0;
        let mut pos = mmap.1;

        loop {
            if pos >= data.len() {
                self.eof = true;
                // EOF: return remaining sequence.
                if self.current_id.is_empty() {
                    return Ok(None);
                }
                if !self.current_sequence.is_empty() {
                    let seq = self.current_fragment()?;
                    self.current_sequence.clear();
                    self.current_id.clear();
                    return Ok(Some(seq));
                }
                return Err(empty_sequence_error(&self.current_id));
            }
            let b = data[pos];
            match b {
                b'>' => {
                    // Hit next sequence: return current accumulation, switch to new ID.
                    let seq = self.current_fragment()?;
                    self.offset = 0;
                    self.current_sequence = String::new();
                    // Read identifier (up to newline, \r-tolerant).
                    pos += 1;
                    let (identifier, next_pos) = parse_mmap_header(&data, pos)?;
                    self.register_id(identifier)?;
                    pos = next_pos;
                    self.mmap = Some((data, pos));
                    return Ok(Some(seq));
                }
                b'\n' | b'\r' => {
                    // Newline (\n or \r, CRLF-tolerant) skipped, not accumulated.
                    pos += 1;
                }
                _ => {
                    self.validate_sequence_byte(b)?;
                    if self.current_sequence.len() == self.fragment_length {
                        // Chop.
                        let seq = self.current_fragment()?;
                        self.offset = self.offset + self.fragment_length - self.overlap;
                        let keep = self.fragment_length - self.overlap;
                        let tail: String = self.current_sequence[keep..].to_string();
                        self.current_sequence = tail;
                        self.current_sequence.push(b as char);
                        self.mmap = Some((data, pos + 1));
                        return Ok(Some(seq));
                    } else {
                        self.current_sequence.push(b as char);
                        pos += 1;
                    }
                }
            }
        }
    }

    /// reader path (gzip or in-memory Cursor).
    fn next_sequence_reader(&mut self) -> std::io::Result<Option<BioSequence>> {
        loop {
            let ch = self.reader.read_byte()?;
            match ch {
                None => {
                    self.eof = true;
                    if self.current_id.is_empty() {
                        return Ok(None);
                    }
                    if !self.current_sequence.is_empty() {
                        let seq = self.current_fragment()?;
                        self.current_sequence.clear();
                        self.current_id.clear();
                        return Ok(Some(seq));
                    }
                    return Err(empty_sequence_error(&self.current_id));
                }
                Some(b'>') => {
                    let seq = self.current_fragment()?;
                    self.offset = 0;
                    self.current_sequence = String::new();
                    let identifier = self.reader.read_header_id()?;
                    self.register_id(identifier)?;
                    return Ok(Some(seq));
                }
                Some(b'\n') | Some(b'\r') => {}
                Some(b) => {
                    self.validate_sequence_byte(b)?;
                    if self.current_sequence.len() == self.fragment_length {
                        let seq = self.current_fragment()?;
                        self.offset = self.offset + self.fragment_length - self.overlap;
                        let keep = self.fragment_length - self.overlap;
                        let tail: String = self.current_sequence[keep..].to_string();
                        self.current_sequence = tail;
                        self.current_sequence.push(b as char);
                        return Ok(Some(seq));
                    } else {
                        self.current_sequence.push(b as char);
                    }
                }
            }
        }
    }

    fn current_fragment(&self) -> std::io::Result<BioSequence> {
        if self.current_sequence.is_empty() {
            return Err(empty_sequence_error(&self.current_id));
        }
        Ok(BioSequence::new(
            format!("{}_{}", self.current_id, self.offset),
            self.current_sequence.to_uppercase(),
        ))
    }

    fn register_id(&mut self, identifier: String) -> std::io::Result<()> {
        if !self.seen_ids.insert(identifier.clone()) {
            return Err(invalid_data(format!(
                "duplicate FASTA identifier '{identifier}'"
            )));
        }
        self.current_id = identifier;
        Ok(())
    }

    fn validate_sequence_byte(&self, byte: u8) -> std::io::Result<()> {
        if is_valid_sequence_byte(byte) {
            return Ok(());
        }
        Err(invalid_data(format!(
            "FASTA record '{}' contains invalid sequence byte 0x{byte:02X}; \
             only ASCII letters, '*' and '-' are accepted",
            self.current_id
        )))
    }
}

fn parse_mmap_header(data: &[u8], mut pos: usize) -> std::io::Result<(String, usize)> {
    let id_start = pos;
    while pos < data.len() && data[pos] != b'\n' && data[pos] != b'\r' {
        pos += 1;
    }
    let header = std::str::from_utf8(&data[id_start..pos]).map_err(|error| {
        invalid_data(format!(
            "FASTA header contains invalid UTF-8 at byte {id_start}: {error}"
        ))
    })?;
    let (identifier, _) = parse_header_fields(header)?;
    while pos < data.len() && (data[pos] == b'\n' || data[pos] == b'\r') {
        pos += 1;
    }
    Ok((identifier, pos))
}

fn empty_sequence_error(identifier: &str) -> std::io::Error {
    invalid_data(format!("FASTA record '{identifier}' has an empty sequence"))
}

/// Internal extension: byte-by-byte read + read identifier.
trait ReadByte {
    fn read_byte(&mut self) -> std::io::Result<Option<u8>>;
    fn read_header_id(&mut self) -> std::io::Result<String>;
}

impl<R: BufRead + ?Sized> ReadByte for R {
    fn read_byte(&mut self) -> std::io::Result<Option<u8>> {
        let mut buf = [0u8; 1];
        match self.read(&mut buf)? {
            0 => Ok(None),
            _ => Ok(Some(buf[0])),
        }
    }

    fn read_header_id(&mut self) -> std::io::Result<String> {
        let mut line = String::new();
        self.read_line(&mut line)?;
        let (identifier, _) = parse_header_fields(line.trim_end())?;
        Ok(identifier)
    }
}

#[cfg(test)]
mod tests {
    use super::SequenceChopper;
    use std::io::Cursor;

    fn chopper(data: &[u8], frag: usize, overlap: usize) -> SequenceChopper {
        let reader: Box<dyn std::io::BufRead> = Box::new(Cursor::new(data.to_vec()));
        SequenceChopper::new(reader, frag, overlap).unwrap()
    }

    #[test]
    fn chop_basic() {
        // 10-base sequence, frag=4, overlap=2.
        let mut c = chopper(b">chr1\nACGTACGTAC\n", 4, 2);
        let s1 = c.next_sequence().unwrap().unwrap();
        assert_eq!(s1.identifier, "chr1_0");
        assert_eq!(s1.sequence, "ACGT");
        let s2 = c.next_sequence().unwrap().unwrap();
        assert_eq!(s2.identifier, "chr1_2"); // offset step = 4-2 = 2
        assert_eq!(s2.sequence, "GTAC");
    }

    #[test]
    fn offset_advances_by_fragment_minus_overlap() {
        let mut c = chopper(b">x\nAAAAAAAAAA\n", 5, 1);
        let s1 = c.next_sequence().unwrap().unwrap();
        assert_eq!(s1.identifier, "x_0");
        let s2 = c.next_sequence().unwrap().unwrap();
        assert_eq!(s2.identifier, "x_4"); // 5-1=4
    }
}
