//! nlr-report — run statistics report (`--stats`).
//!
//! Emits global + per-chromosome + motif dimension statistics, in TSV format (deterministic order).
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::Write;

/// Run statistics result (for TSV output).
pub struct RunStats {
    /// Sequence count.
    pub seq_count: usize,
    /// Total motif count.
    pub motif_total: usize,
    /// Total NLR loci count.
    pub nlr_total: usize,
    /// Complete DNA-mode NLR loci (P-loop plus at least one LRR motif).
    pub nlr_complete: usize,
    /// Per-chromosome statistics: seq -> (motif count, NLR count, complete NLR count).
    pub per_chromosome: BTreeMap<String, (usize, usize, usize)>,
    /// Hit count per motif id.
    pub motif_counts: BTreeMap<u8, usize>,
    /// NLR 结构域类型计数（domain_string -> count）。
    pub nlr_type_counts: BTreeMap<String, usize>,
    /// 生物学家族计数（TNL/RNL/CNL/NB-only/TIR-only/Other）。
    pub nlr_family_counts: BTreeMap<String, usize>,
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
    let nlr_complete = nlrs.iter().filter(|list| list.is_complete_nlr(def)).count();
    let mut per_chromosome: BTreeMap<String, (usize, usize, usize)> = BTreeMap::new();
    for (seq, motifs) in motifs_by_seq {
        let on_seq = nlrs
            .iter()
            .filter(|list| list.sequence_name() == seq.as_str());
        let nlr_for_seq = on_seq.clone().count();
        let complete_for_seq = on_seq.filter(|list| list.is_complete_nlr(def)).count();
        per_chromosome.insert(seq.clone(), (motifs.len(), nlr_for_seq, complete_for_seq));
    }

    let mut motif_counts: BTreeMap<u8, usize> = BTreeMap::new();
    for v in motifs_by_seq.values() {
        for m in v {
            *motif_counts.entry(m.id).or_insert(0) += 1;
        }
    }

    let mut nlr_type_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut nlr_family_counts: BTreeMap<String, usize> = BTreeMap::new();
    for list in nlrs {
        *nlr_type_counts.entry(list.domain_string(def)).or_default() += 1;
        *nlr_family_counts
            .entry(list.family(def).to_string())
            .or_default() += 1;
    }

    RunStats {
        seq_count,
        motif_total,
        nlr_total,
        nlr_complete,
        per_chromosome,
        motif_counts,
        nlr_type_counts,
        nlr_family_counts,
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
    for (family, count) in &stats.nlr_family_counts {
        writeln!(w, "family\t{}\t{}", family, count)?;
    }
    Ok(())
}

/// 单条染色体的统计信息（summary.json 用）。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromosomeStats {
    pub motifs: usize,
    pub nlrs: usize,
    /// Complete DNA-mode NLR count on this sequence.
    #[serde(default)]
    pub complete_nlrs: usize,
    /// 序列长度（bp）；旧版 summary.json 缺失时按 0 处理。
    #[serde(default)]
    pub length: u64,
}

/// 单个 NLR 位点的结构化记录（summary.json 用）。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocusRecord {
    pub name: String,
    pub seq: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub class: String,
    pub family: String,
    /// DNA-mode complete status (P-loop plus at least one LRR motif).
    #[serde(default)]
    pub complete: bool,
    pub has_stop: bool,
    pub motif_count: usize,
}

/// 面向 GUI/Web/R 的机器可读汇总（summary.json）。
///
/// 由 [`build_summary_report`] 从一次运行的 motif 与位点结果直接构造，
/// 字段顺序与值均确定，便于下游做 diff 与解析。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SummaryReport {
    pub sequence_count: usize,
    pub motif_total: usize,
    pub nlr_total: usize,
    /// Complete DNA-mode NLR count.
    #[serde(default)]
    pub nlr_complete: usize,
    pub per_chromosome: BTreeMap<String, ChromosomeStats>,
    pub motif_counts: BTreeMap<u8, usize>,
    pub nlr_type_counts: BTreeMap<String, usize>,
    pub nlr_family_counts: BTreeMap<String, usize>,
    pub loci: Vec<LocusRecord>,
}

/// 从运行结果构造 summary.json 的内容。
pub fn build_summary_report(
    motifs_by_seq: &std::collections::HashMap<String, Vec<Motif>>,
    nlrs: &[MotifList],
    def: &AnnotatorSignatureDefinition,
) -> SummaryReport {
    build_summary_report_with_lengths(motifs_by_seq, nlrs, def, &BTreeMap::new())
}

/// 构造带真实序列长度的 summary.json；用于 GUI 染色体全景图。
pub fn build_summary_report_with_lengths(
    motifs_by_seq: &std::collections::HashMap<String, Vec<Motif>>,
    nlrs: &[MotifList],
    def: &AnnotatorSignatureDefinition,
    sequence_lengths: &BTreeMap<String, u64>,
) -> SummaryReport {
    let stats = collect(motifs_by_seq, nlrs, def);

    let mut nlr_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut complete_counts: BTreeMap<String, usize> = BTreeMap::new();
    for list in nlrs {
        let seq = list.sequence_name().to_string();
        *nlr_counts.entry(seq.clone()).or_default() += 1;
        if list.is_complete_nlr(def) {
            *complete_counts.entry(seq).or_default() += 1;
        }
    }

    let mut per_chromosome: BTreeMap<String, ChromosomeStats> = sequence_lengths
        .iter()
        .map(|(seq, length)| {
            (
                seq.clone(),
                ChromosomeStats {
                    motifs: motifs_by_seq.get(seq).map_or(0, Vec::len),
                    nlrs: nlr_counts.get(seq).copied().unwrap_or(0),
                    complete_nlrs: complete_counts.get(seq).copied().unwrap_or(0),
                    length: *length,
                },
            )
        })
        .collect();

    // 保留没有长度信息的旧调用路径，同时兼容只出现在 motif 表中的序列。
    for (seq, motifs) in motifs_by_seq {
        per_chromosome
            .entry(seq.clone())
            .or_insert_with(|| ChromosomeStats {
                motifs: motifs.len(),
                nlrs: nlr_counts.get(seq).copied().unwrap_or(0),
                complete_nlrs: complete_counts.get(seq).copied().unwrap_or(0),
                length: 0,
            });
    }

    let loci = nlrs
        .iter()
        .map(|list| {
            let (start, end) = list.span();
            LocusRecord {
                name: list.name.clone(),
                seq: list.sequence_name().to_string(),
                start,
                end,
                strand: list.strand().symbol(),
                class: list.domain_string(def),
                family: list.family(def).to_string(),
                complete: list.is_complete_nlr(def),
                has_stop: list.has_stop_codon(),
                motif_count: list.motifs.len(),
            }
        })
        .collect();

    SummaryReport {
        sequence_count: per_chromosome.len(),
        motif_total: stats.motif_total,
        nlr_total: stats.nlr_total,
        nlr_complete: stats.nlr_complete,
        per_chromosome,
        motif_counts: stats.motif_counts,
        nlr_type_counts: stats.nlr_type_counts,
        nlr_family_counts: stats.nlr_family_counts,
        loci,
    }
}

/// 把 summary.json 内容序列化为 JSON 并写入目标 writer。
pub fn write_summary_json<W: Write>(w: &mut W, report: &SummaryReport) -> std::io::Result<()> {
    let bytes = serde_json::to_vec_pretty(report)
        .map_err(|e| std::io::Error::other(format!("serialize summary.json: {e}")))?;
    w.write_all(&bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nlr_core::motif::Motif;
    use nlr_core::strand::Strand;

    fn motif(id: u8, start: u64) -> Motif {
        let mut hit = Motif::new_protein(id, "chr1".to_string(), 1, "A".to_string(), 1e-6);
        hit.set_dna("chr1".to_string(), start, 20_000, 0, Strand::Forward);
        hit
    }

    #[test]
    fn summary_includes_sequences_without_motif_hits() {
        let motifs = std::collections::HashMap::new();
        let def = AnnotatorSignatureDefinition::new();
        let lengths = BTreeMap::from([("chr1".to_string(), 123_456_u64)]);

        let report = build_summary_report_with_lengths(&motifs, &[], &def, &lengths);

        assert_eq!(report.sequence_count, 1);
        assert_eq!(report.per_chromosome["chr1"].length, 123_456);
        assert_eq!(report.per_chromosome["chr1"].motifs, 0);
        assert_eq!(report.per_chromosome["chr1"].nlrs, 0);
    }

    #[test]
    fn chromosome_length_is_optional_for_older_summary_json() {
        let stats: ChromosomeStats = serde_json::from_str(r#"{"motifs": 12, "nlrs": 3}"#).unwrap();
        assert_eq!(stats.length, 0);
        assert_eq!(stats.motifs, 12);
        assert_eq!(stats.nlrs, 3);
        assert_eq!(stats.complete_nlrs, 0);
    }

    #[test]
    fn complete_counts_are_consistent_across_tsv_and_summary() {
        let def = AnnotatorSignatureDefinition::new();
        let complete = MotifList::new("chr1_nlr1".to_string(), vec![motif(1, 100), motif(9, 300)]);
        let partial = MotifList::new("chr1_nlr2".to_string(), vec![motif(6, 500)]);
        let nlrs = vec![complete, partial];
        let motifs = std::collections::HashMap::from([(
            "chr1".to_string(),
            vec![motif(1, 100), motif(9, 300), motif(6, 500)],
        )]);
        let lengths = BTreeMap::from([("chr1".to_string(), 10_000_u64)]);

        let report = build_summary_report_with_lengths(&motifs, &nlrs, &def, &lengths);
        assert_eq!(report.nlr_total, 2);
        assert_eq!(report.nlr_complete, 1);
        assert_eq!(report.per_chromosome["chr1"].complete_nlrs, 1);
        assert!(report.loci[0].complete);
        assert!(!report.loci[1].complete);

        let stats = collect(&motifs, &nlrs, &def);
        let mut tsv = Vec::new();
        write_tsv(&mut tsv, &stats).unwrap();
        let text = String::from_utf8(tsv).unwrap();
        assert!(text.contains("global\tnlr_complete\t1\n"));
        assert!(text.contains("chromosome\tchr1\t3\t2\t1\n"));
    }

    #[test]
    fn old_summary_without_complete_fields_deserializes() {
        let report: SummaryReport = serde_json::from_str(
            r#"{
                "sequence_count": 1,
                "motif_total": 2,
                "nlr_total": 1,
                "per_chromosome": {"chr1": {"motifs": 2, "nlrs": 1, "length": 100}},
                "motif_counts": {"1": 1},
                "nlr_type_counts": {"NBARC": 1},
                "nlr_family_counts": {"NB-only": 1},
                "loci": [{
                    "name": "chr1_nlr1",
                    "seq": "chr1",
                    "start": 1,
                    "end": 2,
                    "strand": "+",
                    "class": "NBARC",
                    "family": "NB-only",
                    "has_stop": false,
                    "motif_count": 1
                }]
            }"#,
        )
        .unwrap();
        assert_eq!(report.nlr_complete, 0);
        assert_eq!(report.per_chromosome["chr1"].complete_nlrs, 0);
        assert!(!report.loci[0].complete);
    }
}
