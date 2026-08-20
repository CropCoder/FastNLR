//! nlr-report — run statistics report (`--stats`).
//!
//! Emits global + per-chromosome + motif dimension statistics, in TSV format (deterministic order).

use std::io::Write;
use std::collections::BTreeMap;
use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;

/// Run statistics result (for TSV output).
pub struct RunStats {
    /// Sequence count.
    pub seq_count: usize,
    /// Total motif count.
    pub motif_total: usize,
    /// Total NLR loci count.
    pub nlr_total: usize,
    /// Complete NLR count.
    pub nlr_complete: usize,
    /// Per-chromosome statistics: seq -> (motif count, NLR count, complete NLR count).
    pub per_chromosome: BTreeMap<String, (usize, usize, usize)>,
    /// Hit count per motif id.
    pub motif_counts: BTreeMap<u8, usize>,
}

/// Aggregate statistics from run results.
pub fn collect(
    motifs_by_seq: &std::collections::HashMap<String, Vec<Motif>>,
    nlrs: &[MotifList],
    def: &AnnotatorSignatureDefinition,
) -> RunStats {
    let seq_count = motifs_by_seq.len();
    let motif_total: usize = motifs_by_seq.values().map(|v| v.len()).sum();
    let nlr_total = nlrs.len();
    let nlr_complete = nlrs.iter().filter(|l| l.is_complete_nlr(def)).count();

    let mut per_chromosome: BTreeMap<String, (usize, usize, usize)> = BTreeMap::new();
    for (seq, motifs) in motifs_by_seq {
        let nlr_for_seq = nlrs.iter().filter(|l| l.sequence_name() == seq.as_str()).count();
        let complete_for_seq = nlrs
            .iter()
            .filter(|l| l.sequence_name() == seq.as_str() && l.is_complete_nlr(def))
            .count();
        per_chromosome.insert(seq.clone(), (motifs.len(), nlr_for_seq, complete_for_seq));
    }

    let mut motif_counts: BTreeMap<u8, usize> = BTreeMap::new();
    for v in motifs_by_seq.values() {
        for m in v {
            *motif_counts.entry(m.id).or_insert(0) += 1;
        }
    }

    RunStats {
        seq_count,
        motif_total,
        nlr_total,
        nlr_complete,
        per_chromosome,
        motif_counts,
    }
}

/// Write the TSV statistics report.
pub fn write_tsv<W: Write>(w: &mut W, stats: &RunStats) -> std::io::Result<()> {
    writeln!(w, "#section\tkey\tvalue")?;
    writeln!(w, "global\tsequence_count\t{}", stats.seq_count)?;
    writeln!(w, "global\tmotif_total\t{}", stats.motif_total)?;
    writeln!(w, "global\tnlr_total\t{}", stats.nlr_total)?;
    writeln!(w, "global\tnlr_complete\t{}", stats.nlr_complete)?;

    for (seq, (motifs, nlrs, complete)) in &stats.per_chromosome {
        writeln!(w, "chromosome\t{}\t{}\t{}\t{}", seq, motifs, nlrs, complete)?;
    }

    for (id, count) in &stats.motif_counts {
        writeln!(w, "motif\tmotif_{}\t{}", id, count)?;
    }
    Ok(())
}
