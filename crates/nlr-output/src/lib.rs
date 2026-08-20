//! nlr-output — output formatting: txt / GFF / BED / motifBED / alignment fasta / loci fasta / TSV.
//!
//! Faithful reimplementation of Java `NLR_Annotator` `write*` methods, with deliberate fixes:
//! 1. `-a` P-loop location logic (`while(!is_ploop)`);
//! 2. `-a` `*`/`X` replacement with `_` (actually applied);
//! 3. `-f` unified clamp on the last contig;
//! 4. `-f` extracts the exact motif span (Java `writeNLRLoci` had an off-by-one that shifted
//!    extraction +1 relative to its own GFF coordinates; this port keeps the extraction
//!    consistent with GFF/BED coordinates — biologically correct);
//! 5. `##date` header set to the current system time (Java was non-deterministic / hardcoded);
//! 6. GFF `##source-version` and source column relabeled `FastNLR`.

use std::io::Write;
use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;

/// Loci report txt (`-o`): no header, 7 columns per NLR.
pub fn write_report_txt<W: Write>(
    w: &mut W,
    nlrs: &[MotifList],
    def: &AnnotatorSignatureDefinition,
) -> std::io::Result<()> {
    for list in nlrs {
        let (start, end) = list.span();
        let strand = list.first_motif().strand.symbol();
        let domains = list.domain_string(def);
        let motif_list = list.motif_list_string();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            list.sequence_name(),
            list.name,
            domains,
            start,
            end,
            strand,
            motif_list
        )?;
    }
    Ok(())
}

/// Loci GFF (`-g`): 6 header lines, start+1 (1-based), end unchanged, feature=NBSLRR.
pub fn write_nlr_gff<W: Write>(
    w: &mut W,
    nlrs: &[MotifList],
    def: &AnnotatorSignatureDefinition,
    date: &str,
    complete_only: bool,
) -> std::io::Result<()> {
    writeln!(w, "##gff-version 2")?;
    writeln!(w, "##source-version FastNLR V1.0")?;
    writeln!(w, "##date {}", date)?;
    writeln!(w, "##Type DNA")?;
    writeln!(w, "#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute")?;

    for list in nlrs {
        if complete_only && !list.is_complete_nlr(def) {
            continue;
        }
        let (start, end) = list.span();
        let strand = list.first_motif().strand.symbol();
        let domains = list.domain_string(def);
        writeln!(
            w,
            "{}\tFastNLR\tNBSLRR\t{}\t{}\t.\t{}\t.\tname={};nlrClass={}",
            list.sequence_name(),
            start + 1, // GFF 1-based
            end,
            strand,
            list.name,
            domains
        )?;
    }
    Ok(())
}

/// Loci BED (`-b`): 12 columns, reverse-strand blocks reversed, colors green/orange/red.
pub fn write_nlr_bed<W: Write>(
    w: &mut W,
    nlrs: &[MotifList],
    def: &AnnotatorSignatureDefinition,
) -> std::io::Result<()> {
    writeln!(w, "#track name=\"NLR_Loci\"")?;
    writeln!(w, "#itemRgb=\"On\"")?;

    for list in nlrs {
        let (start, end) = list.span();
        let strand = list.first_motif().strand.symbol();
        let is_forward = list.is_forward();

        let mut block_starts: Vec<u64> = Vec::with_capacity(list.motifs.len());
        let mut block_sizes: Vec<u64> = Vec::with_capacity(list.motifs.len());
        for m in &list.motifs {
            let rel = m.dna_start - start;
            let size = m.dna_end - m.dna_start;
            if is_forward {
                block_starts.push(rel);
                block_sizes.push(size);
            } else {
                block_starts.insert(0, rel); // reverse strand: reversed order
                block_sizes.insert(0, size);
            }
        }

        let color = if list.has_stop_codon() {
            "255,0,0" // red (stop codon takes priority)
        } else if list.is_complete_nlr(def) {
            "0,255,0" // green (complete)
        } else {
            "255,128,0" // orange (partial)
        };

        let block_count = list.motifs.len();
        let sizes = block_sizes.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
        let starts = block_starts.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");

        writeln!(
            w,
            "{}\t{}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            list.sequence_name(),
            start,
            end,
            list.name,
            strand,
            start,
            end,
            color,
            block_count,
            sizes,
            starts
        )?;
    }
    Ok(())
}

/// Motif BED (`-m`): one row per motif, 9 columns, score column holds pvalue.
pub fn write_motif_bed<W: Write>(
    w: &mut W,
    motifs: &[Motif],
    def: &AnnotatorSignatureDefinition,
    annotate_stop: bool,
) -> std::io::Result<()> {
    writeln!(w, "#track name=\"NLR_Motifs\"")?;
    writeln!(w, "#itemRgb=\"On\"")?;

    for m in motifs {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            m.dna_sequence_id.as_deref().unwrap_or(""),
            m.dna_start,
            m.dna_end,
            nlr_core::signature_def::motif_id_str(m.id),
            nlr_core::motif::format_double_java(m.pvalue),
            m.strand.symbol(),
            m.dna_start,
            m.dna_end,
            def.color_rgb(m.id)
        )?;
        if annotate_stop && m.has_stop() {
            writeln!(
                w,
                "{}\t{}\t{}\tSTOP\t0\t{}\t{}\t{}\t0,0,0",
                m.dna_sequence_id.as_deref().unwrap_or(""),
                m.dna_start,
                m.dna_end,
                m.strand.symbol(),
                m.dna_start,
                m.dna_end
            )?;
        }
    }
    Ok(())
}

/// NB-ARC multiple alignment fasta (`-a`): fixed P-loop location and `*`/`X` replacement.
pub fn write_nbarc_alignment_fasta<W: Write>(
    w: &mut W,
    nlrs: &[MotifList],
    def: &AnnotatorSignatureDefinition,
    include_ced4: bool,
) -> std::io::Result<()> {
    let motif_order = def.nbarc_motif_order();

    if include_ced4 {
        // CED4 reference row: 85-char reference sequence + 80-char dash padding, on one line
        // (matches Java's single literal string for this row).
        writeln!(
            w,
            ">NP_001021202.1\nFLHGRAGSGKSVIASQALSKS-----------------------------TLFVFDDVVQEETIRLRLRCLVTTRDVEISNAASQ{}",
            "-".repeat(80)
        )?;
    }

    for list in nlrs {
        if !list.is_complete_nlr(def) {
            continue;
        }

        // Locate P-loop (fixed Java `while(isPloop)` to `while(!is_ploop)`).
        let ploop_idx = match list.motifs.iter().position(|m| def.is_ploop(m.id)) {
            Some(i) => i,
            None => continue,
        };

        let mut sequence = list.motifs[ploop_idx].protein_sequence.clone();

        // Starting one past ploop, concatenate along motif_order.
        let mut idx = ploop_idx + 1;
        for &expected in motif_order.iter().skip(1) {
            // If the current motif matches expected, append; otherwise pad with '-'.
            let matched = list.motifs.get(idx).filter(|m| m.id == expected);
            match matched {
                Some(m) => {
                    sequence.push_str(&m.protein_sequence);
                    idx += 1;
                }
                None => {
                    let gap_len = def.consensus(expected).len();
                    sequence.push_str(&"-".repeat(gap_len));
                }
            }
        }

        // Fix: actually replace stop codons and unknown aa.
        sequence = sequence.replace('*', "_").replace('X', "_");

        writeln!(w, ">{}\n{}", list.name, sequence)?;
    }
    Ok(())
}

/// Loci sequence fasta for a single contig (`-f`): extract ± flanking from the genome,
/// reverse-complement on the reverse strand, unified clamp.
pub fn write_nlr_loci<W: Write>(
    w: &mut W,
    nlrs: &[MotifList],
    genome: &str,
    genome_id: &str,
    flanking: u64,
) -> std::io::Result<()> {
    let genome_len = genome.len() as u64;
    for list in nlrs {
        let (start, end) = list.span();
        let s = start.saturating_sub(flanking);
        let e = (end + flanking).min(genome_len); // unified clamp
        if s >= e {
            continue;
        }
        write_one_locus(w, list, &genome[s as usize..e as usize], genome_id, s, e)?;
    }
    Ok(())
}

/// Loci sequence fasta across multiple contigs (replicates Java `writeNLRLoci` semantics).
///
/// `contigs` is a list of `(identifier, sequence)` in genome FASTA order. NLRs are grouped
/// by sequence name and extracted from the matching contig. This fixes the original Rust bug
/// that only processed the first contig (`seqs.first()`), which dropped or mis-clamped NLRs
/// on every other chromosome.
pub fn write_nlr_loci_all<W: Write>(
    w: &mut W,
    nlrs: &[MotifList],
    contigs: &[(&str, &str)],
    flanking: u64,
) -> std::io::Result<()> {
    for (genome_id, genome) in contigs {
        let genome_len = genome.len() as u64;
        for list in nlrs.iter().filter(|l| l.sequence_name() == *genome_id) {
            let (start, end) = list.span();
            let s = start.saturating_sub(flanking);
            let e = (end + flanking).min(genome_len);
            if s >= e {
                continue;
            }
            write_one_locus(w, list, &genome[s as usize..e as usize], genome_id, s, e)?;
        }
    }
    Ok(())
}

/// Write a single locus record (shared by single- and multi-contig entry points).
fn write_one_locus<W: Write>(
    w: &mut W,
    list: &MotifList,
    sub: &str,
    genome_id: &str,
    s: u64,
    e: u64,
) -> std::io::Result<()> {
    let seq_str = if list.is_forward() {
        sub.to_string()
    } else {
        // reverse complement
        reverse_complement(sub)
    };
    let strand = list.strand();
    let desc = format!(
        "{} {}-{} strand:{} {}",
        genome_id,
        s,
        e,
        strand.symbol(),
        list.motif_list_string()
    );
    writeln!(w, ">{} {}", list.name, desc)?;
    // wrap at 100 chars per line
    for chunk in seq_str.as_bytes().chunks(100) {
        writeln!(w, "{}", String::from_utf8_lossy(chunk))?;
    }
    Ok(())
}

/// Reverse complement (output-only; inlined here to avoid a cross-crate dependency).
fn reverse_complement(s: &str) -> String {
    let comp: Vec<char> = s
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
            other => other,
        })
        .collect();
    comp.into_iter().rev().collect()
}

/// Motif TSV export (`-c` export): one row per motif, 11 fields.
pub fn export_motifs<W: Write>(w: &mut W, motifs: &[Motif]) -> std::io::Result<()> {
    for m in motifs {
        writeln!(w, "{}", m.export_string())?;
    }
    Ok(())
}

/// Motif TSV import (`-c` import): parse line by line.
pub fn import_motifs(lines: &str) -> Vec<Motif> {
    lines
        .lines()
        .filter_map(Motif::from_export_line)
        .collect()
}
