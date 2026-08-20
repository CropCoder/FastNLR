//! Sequence chopper (streaming slicing with overlap).
//!
//! Reimplementation of Java `SequenceChopper`:
//! - reads character by character, accumulating up to `fragmentLength` then chopping;
//! - fragment naming `{id}_{offset}`, where offset is the 0-based start of the fragment, step `fragmentLength - overlap`;
//! - sequence is uppercased; newlines are skipped; header takes only the first whitespace-delimited token with `>` stripped.
//!
//! Note: Java `readIdentifier` can infinite-loop at EOF (`(char)baseCharacter != '\n'` never terminates);
//! the Rust version adds an EOF guard.

use crate::translate::BioSequence;
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
}

impl SequenceChopper {
    /// Construct from file, auto-sniffing gzip.
    /// Non-gzip files take the mmap path (zero-copy slicing); gzip takes the streaming decompression path.
    pub fn from_file(
        path: &std::path::Path,
        fragment_length: usize,
        overlap: usize,
    ) -> std::io::Result<Self> {
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
            while pos < data.len() && data[pos] != b'\n' && data[pos] != b'\r' {
                current_id.push(data[pos] as char);
                pos += 1;
            }
            // Skip newline (\n or \r\n).
            while pos < data.len() && (data[pos] == b'\n' || data[pos] == b'\r') {
                pos += 1;
            }
            // Take only the first whitespace-delimited token.
            current_id = current_id
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
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
        })
    }

    /// Construct from any BufRead.
    pub fn new(
        mut reader: Box<dyn BufRead>,
        fragment_length: usize,
        overlap: usize,
    ) -> std::io::Result<Self> {
        // Skip leading non-'>' content.
        let mut c = reader.read_byte();
        while let Some(ch) = c {
            if ch == b'>' {
                break;
            }
            c = reader.read_byte();
        }
        // Now c == Some(b'>'); read identifier.
        let current_id = if c == Some(b'>') {
            reader.read_identifier()?
        } else {
            String::new() // empty file
        };

        Ok(SequenceChopper {
            reader,
            fragment_length,
            overlap,
            current_id,
            current_sequence: String::new(),
            offset: 0,
            eof: false,
            mmap: None,
        })
    }

    /// Read the next fragment; returns None at EOF.
    pub fn next_sequence(&mut self) -> Option<BioSequence> {
        if self.eof {
            return None;
        }
        if self.mmap.is_some() {
            return self.next_sequence_mmap();
        }
        self.next_sequence_reader()
    }

    /// mmap path: bytes come from the mmap slice (saves syscalls); accumulation logic matches the reader path.
    fn next_sequence_mmap(&mut self) -> Option<BioSequence> {
        // Take out the current mapping and position.
        let mmap = self.mmap.take()?;
        let data = mmap.0;
        let mut pos = mmap.1;

        loop {
            if pos >= data.len() {
                self.eof = true;
                // EOF: return remaining sequence.
                if !self.current_sequence.is_empty() {
                    let seq = BioSequence::new(
                        format!("{}_{}", self.current_id, self.offset),
                        self.current_sequence.to_uppercase(),
                    );
                    self.current_sequence.clear();
                    return Some(seq);
                }
                return None;
            }
            let b = data[pos];
            match b {
                b'>' => {
                    // Hit next sequence: return current accumulation, switch to new ID.
                    let seq = BioSequence::new(
                        format!("{}_{}", self.current_id, self.offset),
                        self.current_sequence.to_uppercase(),
                    );
                    self.offset = 0;
                    self.current_sequence = String::new();
                    // Read identifier (up to newline, \r-tolerant).
                    pos += 1;
                    let id_start = pos;
                    while pos < data.len() && data[pos] != b'\n' && data[pos] != b'\r' {
                        pos += 1;
                    }
                    let id = String::from_utf8_lossy(&data[id_start..pos])
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_string();
                    self.current_id = id;
                    // Skip newline (\n or \r\n).
                    while pos < data.len() && (data[pos] == b'\n' || data[pos] == b'\r') {
                        pos += 1;
                    }
                    self.mmap = Some((data, pos));
                    return Some(seq);
                }
                b'\n' | b'\r' => {
                    // Newline (\n or \r, CRLF-tolerant) skipped, not accumulated.
                    pos += 1;
                }
                _ => {
                    if self.current_sequence.len() == self.fragment_length {
                        // Chop.
                        let seq = BioSequence::new(
                            format!("{}_{}", self.current_id, self.offset),
                            self.current_sequence.to_uppercase(),
                        );
                        self.offset = self.offset + self.fragment_length - self.overlap;
                        let keep = self.fragment_length - self.overlap;
                        let tail: String = self.current_sequence[keep..].to_string();
                        self.current_sequence = tail;
                        self.current_sequence.push(b as char);
                        self.mmap = Some((data, pos + 1));
                        return Some(seq);
                    } else {
                        self.current_sequence.push(b as char);
                        pos += 1;
                    }
                }
            }
        }
    }

    /// reader path (gzip or in-memory Cursor).
    fn next_sequence_reader(&mut self) -> Option<BioSequence> {
        loop {
            let ch = self.reader.read_byte();
            match ch {
                None => {
                    self.eof = true;
                    if !self.current_sequence.is_empty() {
                        let seq = BioSequence::new(
                            format!("{}_{}", self.current_id, self.offset),
                            self.current_sequence.to_uppercase(),
                        );
                        self.current_sequence.clear();
                        return Some(seq);
                    }
                    return None;
                }
                Some(b'>') => {
                    let seq = BioSequence::new(
                        format!("{}_{}", self.current_id, self.offset),
                        self.current_sequence.to_uppercase(),
                    );
                    self.offset = 0;
                    self.current_sequence = String::new();
                    if let Ok(id) = self.reader.read_identifier() {
                        self.current_id = id;
                    }
                    return Some(seq);
                }
                Some(b'\n') | Some(b'\r') => {}
                Some(b) => {
                    if self.current_sequence.len() == self.fragment_length {
                        let seq = BioSequence::new(
                            format!("{}_{}", self.current_id, self.offset),
                            self.current_sequence.to_uppercase(),
                        );
                        self.offset = self.offset + self.fragment_length - self.overlap;
                        let keep = self.fragment_length - self.overlap;
                        let tail: String = self.current_sequence[keep..].to_string();
                        self.current_sequence = tail;
                        self.current_sequence.push(b as char);
                        return Some(seq);
                    } else {
                        self.current_sequence.push(b as char);
                    }
                }
            }
        }
    }
}

/// Internal extension: byte-by-byte read + read identifier.
trait ReadByte {
    fn read_byte(&mut self) -> Option<u8>;
    fn read_identifier(&mut self) -> std::io::Result<String>;
}

impl<R: BufRead + ?Sized> ReadByte for R {
    fn read_byte(&mut self) -> Option<u8> {
        let mut buf = [0u8; 1];
        match self.read(&mut buf) {
            Ok(0) => None,
            Ok(_) => Some(buf[0]),
            Err(_) => None,
        }
    }

    fn read_identifier(&mut self) -> std::io::Result<String> {
        let mut line = String::new();
        self.read_line(&mut line)?;
        let trimmed = line.trim_end();
        let id = trimmed
            .split_whitespace()
            .next()
            .unwrap_or("")
            .trim_start_matches('>')
            .to_string();
        Ok(id)
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
        let s1 = c.next_sequence().unwrap();
        assert_eq!(s1.identifier, "chr1_0");
        assert_eq!(s1.sequence, "ACGT");
        let s2 = c.next_sequence().unwrap();
        assert_eq!(s2.identifier, "chr1_2"); // offset step = 4-2 = 2
        assert_eq!(s2.sequence, "GTAC");
    }

    #[test]
    fn offset_advances_by_fragment_minus_overlap() {
        let mut c = chopper(b">x\nAAAAAAAAAA\n", 5, 1);
        let s1 = c.next_sequence().unwrap();
        assert_eq!(s1.identifier, "x_0");
        let s2 = c.next_sequence().unwrap();
        assert_eq!(s2.identifier, "x_4"); // 5-1=4
    }
}
