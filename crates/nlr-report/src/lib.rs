//! nlr-report — run statistics report (`--stats`).
//!
//! Emits global + per-chromosome + motif dimension statistics, in TSV format (deterministic order).

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;
use serde::Serialize;
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
    let mut per_chromosome: BTreeMap<String, (usize, usize, usize)> = BTreeMap::new();
    for (seq, motifs) in motifs_by_seq {
        let nlr_for_seq = nlrs
            .iter()
            .filter(|l| l.sequence_name() == seq.as_str())
            .count();
        // 完整性判定暂时禁用：第三项 complete 计数固定为 0 占位。
        per_chromosome.insert(seq.clone(), (motifs.len(), nlr_for_seq, 0));
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
    for (seq, (motifs, nlrs, _complete)) in &stats.per_chromosome {
        writeln!(w, "chromosome\t{}\t{}\t{}", seq, motifs, nlrs)?;
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
#[derive(Debug, Clone, Serialize)]
pub struct ChromosomeStats {
    pub motifs: usize,
    pub nlrs: usize,
}

/// 单个 NLR 位点的结构化记录（summary.json 用）。
#[derive(Debug, Clone, Serialize)]
pub struct LocusRecord {
    pub name: String,
    pub seq: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub class: String,
    pub family: String,
    pub has_stop: bool,
    pub motif_count: usize,
}

/// 面向 GUI/Web/R 的机器可读汇总（summary.json）。
///
/// 由 [`build_summary_report`] 从一次运行的 motif 与位点结果直接构造，
/// 字段顺序与值均确定，便于下游做 diff 与解析。
#[derive(Debug, Clone, Serialize)]
pub struct SummaryReport {
    pub sequence_count: usize,
    pub motif_total: usize,
    pub nlr_total: usize,
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
    let stats = collect(motifs_by_seq, nlrs, def);

    let per_chromosome = stats
        .per_chromosome
        .iter()
        .map(|(seq, (motifs, nlrs, _complete))| {
            (
                seq.clone(),
                ChromosomeStats {
                    motifs: *motifs,
                    nlrs: *nlrs,
                },
            )
        })
        .collect();

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
                has_stop: list.has_stop_codon(),
                motif_count: list.motifs.len(),
            }
        })
        .collect();

    SummaryReport {
        sequence_count: stats.seq_count,
        motif_total: stats.motif_total,
        nlr_total: stats.nlr_total,
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
