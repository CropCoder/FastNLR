//! Strand enum.

/// Strand: forward strand / reverse strand.
///
/// Corresponds to Java `Motif.forwardStrand`.
///
/// Note: on both strands, `Motif.dna_start` is always leftmost and `dna_end` always rightmost
/// (i.e. `dna_start < dna_end` always holds). The only special behavior of the reverse strand
/// is that it is **sorted descending** (genomic coordinates decrease when reading in the
/// translation direction N->C); the coordinates themselves are not inverted.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    /// Construct from a boolean.
    #[inline]
    pub fn from_bool(forward: bool) -> Self {
        if forward {
            Strand::Forward
        } else {
            Strand::Reverse
        }
    }

    /// Convert to boolean (true = forward strand).
    #[inline]
    pub fn is_forward(self) -> bool {
        matches!(self, Strand::Forward)
    }

    /// Convert to output symbol `+` / `-`.
    #[inline]
    pub fn symbol(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}
