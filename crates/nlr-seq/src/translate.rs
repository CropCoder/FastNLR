//! Biological sequence with six-frame translation and reverse complement.
//!
//! Reimplementation of Java `BioSequence`:
//! - `translate2protein`: 6 frames (forward strand 0/1/2 + reverse strand 0/1/2), remainder -> end offset table;
//! - naming `{id}_frame+0/1/2`, `{id}_frame-0/1/2`;
//! - reverse complement: case-sensitive, non-standard bases preserved as-is.

use crate::codon::translate_triplet;

/// Biological sequence (identifier + description + sequence string).
///
/// `identifier` is the clean ID (first whitespace-delimited token of the FASTA header with `>` stripped);
/// `description` is the remainder of the header (may be empty), only appended when emitting FASTA, not used in translation naming.
#[derive(Debug, Clone, PartialEq)]
pub struct BioSequence {
    pub identifier: String,
    pub description: String,
    pub sequence: String,
}

impl BioSequence {
    pub fn new(identifier: impl Into<String>, sequence: impl Into<String>) -> Self {
        BioSequence {
            identifier: identifier.into(),
            description: String::new(),
            sequence: sequence.into(),
        }
    }

    /// Sequence length (in bytes).
    #[inline]
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Reverse complement (equivalent to `getReverseComplementarySequence`).
    /// Complements character by character (case-sensitive), then reverses the whole string.
    pub fn reverse_complement(&self) -> String {
        let comp: Vec<char> = self
            .sequence
            .chars()
            .map(|b| match b {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                'a' => 't',
                't' => 'a',
                'g' => 'c',
                'c' => 'g',
                other => other, // non-standard bases preserved as-is
            })
            .collect();
        comp.into_iter().rev().collect()
    }

    /// Six-frame translation (equivalent to `translate2Protein`).
    ///
    /// Returns 6 BioSequence instances, named `{id}_frame+0..+2`, `{id}_frame-0..-2`.
    pub fn translate2protein(&self) -> Vec<BioSequence> {
        let seq = self.sequence.as_bytes();
        let n = seq.len();
        let terminus = n % 3;

        // Remainder -> end offset table (shared by forward and reverse strands).
        // Java: frame0 t0=-2,t1=-3,t2=-4; frame1 t0=-4,t1=-2,t2=-3; frame2 t0=-3,t1=-4,t2=-2
        let end_offset = |frame: usize| -> isize {
            match (terminus, frame) {
                (0, 0) | (1, 1) | (2, 2) => -2,
                (1, 0) | (2, 1) | (0, 2) => -3,
                (2, 0) | (0, 1) | (1, 2) => -4,
                _ => unreachable!(),
            }
        };

        let translate_frame = |bytes: &[u8], frame: usize| -> String {
            let end = end_offset(frame);
            let limit = (n as isize + end).max(0) as usize;
            let mut out = String::with_capacity(limit / 3 + 1);
            let mut i = frame;
            // Java: for(i=frame; i < n+end; i+=3) substring(i,i+3)
            while i < limit && i + 3 <= n {
                out.push(translate_triplet(&bytes[i..i + 3]));
                i += 3;
            }
            out
        };

        // Forward strand: 3 frames.
        let mut frames = Vec::with_capacity(6);
        for frame in 0..3 {
            let aa = translate_frame(seq, frame);
            frames.push(BioSequence::new(
                format!("{}_frame+{}", self.identifier, frame),
                aa,
            ));
        }

        // Reverse strand: 3 frames (same 0/1/2 offset applied to the reverse complement).
        let rc = self.reverse_complement();
        let rc_bytes = rc.as_bytes();
        let rc_n = rc_bytes.len();
        let rc_translate = |bytes: &[u8], frame: usize| -> String {
            let end = end_offset(frame);
            let limit = (rc_n as isize + end).max(0) as usize;
            let mut out = String::with_capacity(limit / 3 + 1);
            let mut i = frame;
            while i < limit && i + 3 <= rc_n {
                out.push(translate_triplet(&bytes[i..i + 3]));
                i += 3;
            }
            out
        };
        for frame in 0..3 {
            let aa = rc_translate(rc_bytes, frame);
            frames.push(BioSequence::new(
                format!("{}_frame-{}", self.identifier, frame),
                aa,
            ));
        }

        frames
    }

    /// Whether this is DNA (equivalent to Java `isDNA`: ATGCN ratio in the first 500 characters > 50%).
    pub fn is_dna(&self) -> bool {
        let sample = &self.sequence[..self.sequence.len().min(500)];
        if sample.is_empty() {
            return false;
        }
        let count = sample
            .bytes()
            .filter(|b| matches!(b, b'A' | b'T' | b'G' | b'C' | b'N' | b'a' | b't' | b'g' | b'c' | b'n'))
            .count();
        count as f64 / sample.len() as f64 > 0.5
    }

    /// FASTA string (single-line sequence, equivalent to `getFastaString`).
    pub fn fasta_string(&self) -> String {
        if self.description.is_empty() {
            format!(">{}\n{}\n", self.identifier, self.sequence)
        } else {
            format!(">{} {}\n{}\n", self.identifier, self.description, self.sequence)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::BioSequence;

    #[test]
    fn reverse_complement_preserves_non_standard() {
        let s = BioSequence::new("x", "ATGCNn");
        // A->T, T->A, G->C, C->G, N->N, n->n, then reversed.
        assert_eq!(s.reverse_complement(), "nNGCAT");
    }

    #[test]
    fn six_frame_naming() {
        let s = BioSequence::new("chr1_0", "ATGAAATTT");
        let frames = s.translate2protein();
        assert_eq!(frames.len(), 6);
        assert_eq!(frames[0].identifier, "chr1_0_frame+0");
        assert_eq!(frames[3].identifier, "chr1_0_frame-0");
    }

    #[test]
    fn six_frame_length_handling() {
        // Length 8 (%3=2); frame0 should translate floor(8/3)=2 complete codons.
        let s = BioSequence::new("x", "ATGAAATT");
        let frames = s.translate2protein();
        // frame+0: ATG AAA -> MK (2 codons)
        assert_eq!(frames[0].sequence, "MK");
        // frame+1: TGA AAT -> *N (2 codons, starting at index 1)
        assert_eq!(frames[1].sequence, "*N");
        // frame+2: GAA ATT -> EI (2 codons, starting at index 2)
        assert_eq!(frames[2].sequence, "EI");
    }
}
