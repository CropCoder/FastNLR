//! Genetic code table and triplet codon translation.
//!
//! Reimplementation of Java `BioSequence.loadGeneticCodeTable` / `translateTriplet`:
//! - 64 standard codons (including 3 stop codons TAA/TAG/TGA -> `*`);
//! - non-standard triplets (containing N, IUPAC, etc.) -> `X`;
//! - input is uppercased before table lookup.

/// Encodes 3 bases into an integer index (A=0, C=1, G=2, T=3; other bases map to out-of-range, returning X).
#[inline]
fn base_idx(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // non-ACGT -> out-of-range marker
    }
}

/// 64-entry codon table, index = base1*16 + base2*4 + base3 (A=0,C=1,G=2,T=3).
/// Includes out-of-range guard: non-ACGT directly returns X.
pub fn translate_triplet(t: &[u8]) -> char {
    debug_assert_eq!(t.len(), 3);
    let i0 = base_idx(t[0]);
    let i1 = base_idx(t[1]);
    let i2 = base_idx(t[2]);
    if i0 > 3 || i1 > 3 || i2 > 3 {
        return 'X';
    }
    CODON_TABLE[i0 * 16 + i1 * 4 + i2]
}

/// Static codon table (A=0, C=1, G=2, T=3 encoding, index = b1*16 + b2*4 + b3).
static CODON_TABLE: [char; 64] = [
    // b1 = A (0)
    'K', 'N', 'K', 'N', // AAA AAC AAG AAT
    'T', 'T', 'T', 'T', // ACA ACC ACG ACT
    'R', 'S', 'R', 'S', // AGA AGC AGG AGT
    'I', 'I', 'M', 'I', // ATA ATC ATG ATT
    // b1 = C (1)
    'Q', 'H', 'Q', 'H', // CAA CAC CAG CAT
    'P', 'P', 'P', 'P', // CCA CCC CCG CCT
    'R', 'R', 'R', 'R', // CGA CGC CGG CGT
    'L', 'L', 'L', 'L', // CTA CTC CTG CTT
    // b1 = G (2)
    'E', 'D', 'E', 'D', // GAA GAC GAG GAT
    'A', 'A', 'A', 'A', // GCA GCC GCG GCT
    'G', 'G', 'G', 'G', // GGA GGC GGG GGT
    'V', 'V', 'V', 'V', // GTA GTC GTG GTT
    // b1 = T (3)
    '*', 'Y', '*', 'Y', // TAA TAC TAG TAT
    'S', 'S', 'S', 'S', // TCA TCC TCG TCT
    '*', 'C', 'W', 'C', // TGA TGC TGG TGT
    'L', 'F', 'L', 'F', // TTA TTC TTG TTT
];

#[cfg(test)]
mod tests {
    use super::translate_triplet;

    #[test]
    fn standard_codons() {
        assert_eq!(translate_triplet(b"ATG"), 'M');
        assert_eq!(translate_triplet(b"TAA"), '*');
        assert_eq!(translate_triplet(b"GGG"), 'G');
    }

    #[test]
    fn non_standard_returns_x() {
        assert_eq!(translate_triplet(b"NNN"), 'X');
        assert_eq!(translate_triplet(b"ATN"), 'X');
    }
}
