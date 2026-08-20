//! Single motif hit type.
//!
//! Faithfully replicates the Java `Motif`:
//! - Holds both protein coordinates (position, protein_sequence) and DNA coordinates (dna_start, dna_end, strand, frame);
//! - `set_dna` maps protein coordinates to genomic coordinates (including reverse strand mirroring);
//! - The sorting rule is asymmetric: forward strand sorted ascending by dna_start, reverse strand descending.

use std::cmp::Ordering;

use crate::signature_def::MotifId;
use crate::strand::Strand;

/// Format an `f64` close to Java's `Double.toString`: shortest round-trippable decimal,
/// plain notation when `1e-3 <= |d| < 1e7`, otherwise scientific notation with uppercase `E`.
///
/// Rust's `{}` / `{:e}` already produce the shortest digit sequence (same as Java); this only
/// picks the notation and uppercases the exponent marker so output matches the original tool
/// (e.g. p-values render as `8.778909458322434E-12`, not `0.000000000008778909458322434`).
pub fn format_double_java(d: f64) -> String {
    if d == 0.0 {
        return "0.0".to_string();
    }
    if d.is_nan() {
        return "NaN".to_string();
    }
    if d.is_infinite() {
        return if d < 0.0 {
            "-Infinity".to_string()
        } else {
            "Infinity".to_string()
        };
    }
    let abs = d.abs();
    if abs >= 1e-3 && abs < 1e7 {
        format!("{}", d)
    } else {
        format!("{:e}", d).replace('e', "E")
    }
}

/// A single motif hit.
#[derive(Debug, Clone, PartialEq)]
pub struct Motif {
    /// Motif id (1..=20).
    pub id: MotifId,
    /// Protein sequence id (translated fragment name, e.g. `chr1_0_frame+0`).
    pub protein_sequence_id: String,
    /// Protein coordinate (1-based semantics; Java stores `q+1` on construction).
    pub position: u64,
    /// The actual matched amino acid sequence.
    pub protein_sequence: String,
    /// p-value.
    pub pvalue: f64,
    /// Score (always 0.0 on the Java scan path; only read from TSV on import).
    pub score: f64,
    /// DNA sequence id (None when unset).
    pub dna_sequence_id: Option<String>,
    /// DNA start (0-based, leftmost; always < dna_end on both strands).
    pub dna_start: u64,
    /// DNA end (0-based, rightmost).
    pub dna_end: u64,
    /// Strand.
    pub strand: Strand,
    /// Reading frame 0/1/2.
    pub frame: u8,
    /// Whether DNA parameters have been set.
    pub dna_parameters_set: bool,
}

impl Motif {
    /// Construct a motif with only protein-side info (equivalent to the Java scan-stage construction).
    pub fn new_protein(
        id: MotifId,
        protein_sequence_id: String,
        position: u64,
        protein_sequence: String,
        pvalue: f64,
    ) -> Self {
        Motif {
            id,
            protein_sequence_id,
            position,
            protein_sequence,
            pvalue,
            score: 0.0,
            dna_sequence_id: None,
            dna_start: 0,
            dna_end: 0,
            strand: Strand::Forward,
            frame: 0,
            dna_parameters_set: false,
        }
    }

    /// Set DNA coordinates (equivalent to Java `Motif.setDNA`).
    ///
    /// - `offset`: 0-based start of the fragment on the chromosome;
    /// - `fragment_length`: fragment length (needed for reverse-strand computation);
    /// - `frame`: reading frame 0/1/2;
    /// - `strand`: strand.
    ///
    /// Forward strand: `dna_start = (position-1)*3 + frame + offset`
    /// Reverse strand: `dna_start = offset + frag_len - ((position+len-1)*3 + frame)`
    ///
    /// On both strands, `dna_start` is always leftmost and `dna_end` always rightmost (`dna_start < dna_end`).
    pub fn set_dna(
        &mut self,
        dna_sequence_id: String,
        offset: u64,
        fragment_length: u64,
        frame: u8,
        strand: Strand,
    ) {
        // Java position is 1-based; the protein sequence length equals the motif length.
        let len = self.protein_sequence.len() as u64;
        let pos_minus_1 = self.position.saturating_sub(1);
        let (start, end) = match strand {
            Strand::Forward => {
                let s = pos_minus_1 * 3 + frame as u64 + offset;
                let e = (self.position + len - 1) * 3 + frame as u64 + offset;
                (s, e)
            }
            Strand::Reverse => {
                let s = offset + fragment_length - ((self.position + len - 1) * 3 + frame as u64);
                let e = offset + fragment_length - (pos_minus_1 * 3 + frame as u64);
                (s, e)
            }
        };
        self.dna_sequence_id = Some(dna_sequence_id);
        self.dna_start = start;
        self.dna_end = end;
        self.strand = strand;
        self.frame = frame;
        self.dna_parameters_set = true;
    }

    /// Whether the sequence contains a stop codon `*`.
    #[inline]
    pub fn has_stop(&self) -> bool {
        self.protein_sequence.contains('*')
    }

    /// Export a TSV string (equivalent to Java `getExportString`, 11 fields).
    pub fn export_string(&self) -> String {
        let mut s = format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            crate::signature_def::motif_id_str(self.id),
            self.protein_sequence_id,
            self.position,
            self.protein_sequence,
            format_double_java(self.pvalue),
            format_double_java(self.score)
        );
        if self.dna_parameters_set {
            s.push_str(&format!(
                "\t{}\t{}\t{}\t{}\t{}",
                self.dna_sequence_id.as_deref().unwrap_or(""),
                self.dna_start,
                self.dna_end,
                self.strand.symbol(),
                self.frame
            ));
        } else {
            s.push_str("\t\t\t\t\t");
        }
        s
    }

    /// Deserialize from a TSV line (equivalent to Java `Motif(String)`, used for -c import).
    pub fn from_export_line(line: &str) -> Option<Self> {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
            return None;
        }
        // Parse the "motif_N" id.
        let id = cols[0]
            .trim_start_matches("motif_")
            .parse::<MotifId>()
            .ok()?;
        let position: u64 = cols[2].parse().ok()?;
        let pvalue: f64 = cols[4].parse().ok()?;
        let score: f64 = cols[5].parse().ok()?;
        let mut m = Motif::new_protein(
            id,
            cols[1].to_string(),
            position,
            cols[3].to_string(),
            pvalue,
        );
        m.score = score;
        // If the last 5 columns exist and are non-empty, restore DNA parameters.
        if cols.len() >= 11 && !cols[6].is_empty() {
            m.dna_sequence_id = Some(cols[6].to_string());
            m.dna_start = cols[7].parse().ok()?;
            m.dna_end = cols[8].parse().ok()?;
            m.strand = if cols[9] == "-" {
                Strand::Reverse
            } else {
                Strand::Forward
            };
            m.frame = cols[10].parse().ok()?;
            m.dna_parameters_set = true;
        }
        Some(m)
    }
}

/// Sort key: only participates in coordinate sorting when DNA parameters are set.
/// (The actual sorting logic is in the Ord impl below; this comment documents the intent.)
impl PartialOrd for Motif {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for Motif {}

impl Ord for Motif {
    /// Replicates the asymmetric sorting of Java `Motif.compareTo`.
    ///
    /// Key point: Java uses `dna_sequence_id` (chromosome id) to determine the same group;
    /// only within the same group does it sort by strand + dnaStart.
    /// Different groups are sorted by `protein_sequence_id`. Do not use protein_sequence_id
    /// to determine the same group (fragment ids are almost unique per motif, which would
    /// prevent same-chromosome motifs from ever being sorted by coordinate).
    fn cmp(&self, other: &Self) -> Ordering {
        if self.dna_parameters_set && other.dna_parameters_set {
            let self_id = self.dna_sequence_id.as_deref().unwrap_or("");
            let other_id = other.dna_sequence_id.as_deref().unwrap_or("");
            if self_id.eq_ignore_ascii_case(other_id) {
                // Same group: first by strand (forward strand first), then by coordinate.
                match (self.strand, other.strand) {
                    (Strand::Forward, Strand::Reverse) => Ordering::Less,
                    (Strand::Reverse, Strand::Forward) => Ordering::Greater,
                    _ => match self.strand {
                        Strand::Forward => self.dna_start.cmp(&other.dna_start),
                        Strand::Reverse => other.dna_start.cmp(&self.dna_start),
                    },
                }
            } else {
                // Different groups: sort by fragment id (replicates Java behavior).
                self.protein_sequence_id.cmp(&other.protein_sequence_id)
            }
        } else {
            // DNA unset: sort by protein id + position.
            if self.protein_sequence_id.eq_ignore_ascii_case(&other.protein_sequence_id) {
                self.position.cmp(&other.position)
            } else {
                self.protein_sequence_id.cmp(&other.protein_sequence_id)
            }
        }
    }
}
