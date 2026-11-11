//! 流式 GFF3 解析与 CDS 转录本模型。
//!
//! 本模块只处理 GFF3 的 `gene`、`mRNA`/`transcript` 和 `CDS` 记录，保留基因名称、
//! biotype、链方向和 CDS phase，并为每个基因选择 CDS 总长度最长的转录本。
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// GFF3 链方向。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn from_symbol(value: &str) -> io::Result<Self> {
        match value {
            "+" => Ok(Strand::Plus),
            "-" => Ok(Strand::Minus),
            other => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid GFF3 strand '{other}'"),
            )),
        }
    }

    pub fn symbol(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
        }
    }
}

/// 一个 CDS 区间，内部坐标为 0-based、half-open。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdsSegment {
    pub seqid: String,
    pub start0: u64,
    pub end0: u64,
    pub phase: Option<u8>,
}

impl CdsSegment {
    pub fn len(&self) -> u64 {
        self.end0.saturating_sub(self.start0)
    }

    pub fn is_empty(&self) -> bool {
        self.start0 >= self.end0
    }
}

/// 基因层注释。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneAnnotation {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub gene_biotype: Option<String>,
    pub seqid: String,
    pub start0: u64,
    pub end0: u64,
    pub strand: Strand,
}

/// 被选中的最长转录本及其 CDS 模型。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptModel {
    pub gene: GeneAnnotation,
    pub transcript_id: String,
    pub transcript_name: Option<String>,
    pub seqid: String,
    pub start0: u64,
    pub end0: u64,
    pub strand: Strand,
    /// 已按转录方向排列的 CDS 片段。
    pub cds: Vec<CdsSegment>,
    pub cds_length_bp: u64,
    pub phase_consistent: bool,
}

impl TranscriptModel {
    /// Bases to skip at the 5' end before translation.
    ///
    /// GFF3 phase on the first CDS in transcript order describes a partial
    /// leading codon. The raw CDS sequence remains unchanged; callers apply
    /// this offset only when translating and mapping protein coordinates.
    pub fn translation_offset(&self) -> usize {
        self.cds
            .first()
            .and_then(|segment| segment.phase)
            .map_or(0, usize::from)
    }

    /// 从按位点方向排列的 contig 序列中拼接 CDS；负链片段逐段反向互补。
    pub fn extract_cds(&self, contig: &str) -> io::Result<String> {
        let mut sequence = String::with_capacity(self.cds_length_bp as usize);
        for segment in &self.cds {
            let start = usize::try_from(segment.start0).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "CDS start exceeds usize")
            })?;
            let end = usize::try_from(segment.end0)
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "CDS end exceeds usize"))?;
            if end > contig.len() || start >= end {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "CDS {}:{}-{} exceeds contig length {}",
                        segment.seqid,
                        segment.start0 + 1,
                        segment.end0,
                        contig.len()
                    ),
                ));
            }
            let slice = &contig[start..end];
            match self.strand {
                Strand::Plus => sequence.push_str(slice),
                Strand::Minus => sequence.push_str(&reverse_complement(slice)),
            }
        }
        Ok(sequence.to_ascii_uppercase())
    }
}

/// 无法进入扫描流程的基因或转录本。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SkippedAnnotation {
    pub gene_id: String,
    pub transcript_id: String,
    pub seqid: String,
    pub reason: String,
}

/// GFF3 解析后的最长转录本集合。
#[derive(Debug, Clone, Default)]
pub struct AnnotationSet {
    pub transcripts: Vec<TranscriptModel>,
    pub skipped: Vec<SkippedAnnotation>,
    /// 至少包含一个 CDS 片段的基因数（选择最长转录本之前）。
    pub genes_with_cds: usize,
}

#[derive(Debug, Clone)]
struct GeneBuilder {
    gene_id: String,
    gene_name: Option<String>,
    gene_biotype: Option<String>,
    seqid: String,
    start0: u64,
    end0: u64,
    strand: Strand,
}

#[derive(Debug, Clone)]
struct TranscriptBuilder {
    transcript_id: String,
    transcript_name: Option<String>,
    gene_id: String,
    seqid: Option<String>,
    start0: Option<u64>,
    end0: Option<u64>,
    strand: Option<Strand>,
    cds: Vec<CdsSegment>,
    error: Option<String>,
}

impl TranscriptBuilder {
    fn placeholder(transcript_id: &str) -> Self {
        TranscriptBuilder {
            transcript_id: transcript_id.to_string(),
            transcript_name: None,
            gene_id: transcript_id.to_string(),
            seqid: None,
            start0: None,
            end0: None,
            strand: None,
            cds: Vec::new(),
            error: None,
        }
    }

    fn observe_cds(&mut self, segment: &CdsSegment, strand: Strand) {
        if self.error.is_some() {
            return;
        }
        if let Some(seqid) = &self.seqid {
            if seqid != &segment.seqid {
                self.error = Some(format!(
                    "transcript {} has CDS on multiple sequences: {} and {}",
                    self.transcript_id, seqid, segment.seqid
                ));
                return;
            }
        } else {
            self.seqid = Some(segment.seqid.clone());
        }
        if let Some(existing) = self.strand {
            if existing != strand {
                self.error = Some(format!(
                    "transcript {} has CDS on conflicting strands",
                    self.transcript_id
                ));
            }
        } else {
            self.strand = Some(strand);
        }
    }
}

/// 读取 GFF3，并返回每个含 CDS 基因的最长转录本。
pub fn read_gff3(path: &Path) -> io::Result<AnnotationSet> {
    let mut file = std::fs::File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    let reader: Box<dyn BufRead> = if n == 2 && magic == [0x1f, 0x8b] {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    read_gff3_from(reader)
}

fn read_gff3_from<R: BufRead>(mut reader: R) -> io::Result<AnnotationSet> {
    let mut genes: HashMap<String, GeneBuilder> = HashMap::new();
    let mut transcripts: HashMap<String, TranscriptBuilder> = HashMap::new();
    let mut line_number = 0usize;
    let mut line = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        line_number += 1;
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == "##FASTA" {
            break;
        }
        if trimmed.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() != 9 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GFF3 line {line_number}: expected 9 columns, found {}",
                    fields.len()
                ),
            ));
        }
        let feature_type = fields[2];
        if !matches!(feature_type, "gene" | "mRNA" | "transcript" | "CDS") {
            continue;
        }
        let seqid = fields[0].to_string();
        let start0 = parse_coordinate(fields[3], line_number)?
            .checked_sub(1)
            .ok_or_else(|| invalid_line(line_number, "start must be >= 1"))?;
        let end0 = parse_coordinate(fields[4], line_number)?;
        if end0 <= start0 {
            return Err(invalid_line(line_number, "end must be greater than start"));
        }
        let strand = match fields[6] {
            "+" | "-" => Strand::from_symbol(fields[6]).map_err(|error| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("GFF3 line {line_number}: {error}"),
                )
            })?,
            "." | "?" => continue,
            other => {
                return Err(invalid_line(
                    line_number,
                    &format!("invalid GFF3 strand '{other}'"),
                ))
            }
        };
        let attributes = parse_attributes(fields[8]);

        match feature_type {
            "gene" => {
                let Some(gene_id) = attributes.get("ID").and_then(|v| v.first()).cloned() else {
                    return Err(invalid_line(line_number, "gene record lacks ID"));
                };
                let mut builder = GeneBuilder {
                    gene_id: gene_id.clone(),
                    gene_name: first_attribute(&attributes, &["Name"]),
                    gene_biotype: first_attribute(
                        &attributes,
                        &["gene_biotype", "biotype", "gene_type"],
                    ),
                    seqid,
                    start0,
                    end0,
                    strand,
                };
                if let Some(previous) = genes.get_mut(&gene_id) {
                    if previous.seqid != builder.seqid || previous.strand != builder.strand {
                        // Multi-part organelle genes can have segments on opposite strands.
                        // The corresponding transcript cannot be represented as one CDS model.
                        continue;
                    }
                    // NCBI GFF3 can represent one gene as several ordered gene segments with
                    // the same public ID. Merge their span rather than rejecting the annotation.
                    previous.start0 = previous.start0.min(builder.start0);
                    previous.end0 = previous.end0.max(builder.end0);
                    if previous.gene_name.is_none() {
                        previous.gene_name = builder.gene_name.take();
                    }
                    if previous.gene_biotype.is_none() {
                        previous.gene_biotype = builder.gene_biotype.take();
                    }
                } else {
                    genes.insert(gene_id, builder);
                }
            }
            "mRNA" | "transcript" => {
                let Some(transcript_id) = attributes.get("ID").and_then(|v| v.first()).cloned()
                else {
                    return Err(invalid_line(
                        line_number,
                        &format!("{feature_type} record lacks ID"),
                    ));
                };
                let gene_id = attributes
                    .get("Parent")
                    .and_then(|v| v.first())
                    .cloned()
                    .unwrap_or_else(|| transcript_id.clone());
                let entry = transcripts
                    .entry(transcript_id.clone())
                    .or_insert_with(|| TranscriptBuilder::placeholder(&transcript_id));
                if entry.seqid.as_ref().is_some_and(|value| value != &seqid)
                    || entry.strand.is_some_and(|value| value != strand)
                {
                    return Err(invalid_line(
                        line_number,
                        &format!("conflicting duplicate transcript ID {transcript_id}"),
                    ));
                }
                entry.transcript_name = first_attribute(&attributes, &["Name"]);
                entry.gene_id = gene_id;
                entry.seqid = Some(seqid);
                entry.start0 = Some(start0);
                entry.end0 = Some(end0);
                entry.strand = Some(strand);
            }
            "CDS" => {
                let Some(parents) = attributes.get("Parent") else {
                    return Err(invalid_line(line_number, "CDS record lacks Parent"));
                };
                let phase = match fields[7] {
                    "." => None,
                    "0" => Some(0),
                    "1" => Some(1),
                    "2" => Some(2),
                    other => {
                        return Err(invalid_line(
                            line_number,
                            &format!("invalid CDS phase '{other}'"),
                        ))
                    }
                };
                for parent in parents {
                    let entry = transcripts
                        .entry(parent.clone())
                        .or_insert_with(|| TranscriptBuilder::placeholder(parent));
                    entry.observe_cds(
                        &CdsSegment {
                            seqid: seqid.clone(),
                            start0,
                            end0,
                            phase,
                        },
                        strand,
                    );
                    entry.cds.push(CdsSegment {
                        seqid: seqid.clone(),
                        start0,
                        end0,
                        phase,
                    });
                }
            }
            _ => unreachable!(),
        }
    }

    finalize_annotation(genes, transcripts)
}

fn finalize_annotation(
    genes: HashMap<String, GeneBuilder>,
    transcripts: HashMap<String, TranscriptBuilder>,
) -> io::Result<AnnotationSet> {
    let mut selected: BTreeMap<String, TranscriptModel> = BTreeMap::new();
    let mut skipped = Vec::new();
    let mut cds_genes = HashSet::new();

    for builder in transcripts.into_values() {
        if builder.cds.is_empty() {
            continue;
        }
        cds_genes.insert(builder.gene_id.clone());
        let result = build_transcript(builder, &genes);
        match result {
            Ok(model) => {
                let replace = selected.get(&model.gene.gene_id).is_none_or(|current| {
                    model.cds_length_bp > current.cds_length_bp
                        || (model.cds_length_bp == current.cds_length_bp
                            && model.transcript_id < current.transcript_id)
                });
                if replace {
                    selected.insert(model.gene.gene_id.clone(), model);
                }
            }
            Err((model, reason)) => skipped.push(SkippedAnnotation {
                gene_id: model.gene_id,
                transcript_id: model.transcript_id,
                seqid: model.seqid,
                reason,
            }),
        }
    }

    if selected.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GFF3 contains no usable gene/transcript/CDS models",
        ));
    }

    let mut transcripts: Vec<TranscriptModel> = selected.into_values().collect();
    transcripts.sort_by(|a, b| {
        (&a.seqid, a.start0, &a.gene.gene_id).cmp(&(&b.seqid, b.start0, &b.gene.gene_id))
    });
    skipped.sort_by(|a, b| {
        (&a.seqid, &a.gene_id, &a.transcript_id).cmp(&(&b.seqid, &b.gene_id, &b.transcript_id))
    });
    Ok(AnnotationSet {
        transcripts,
        skipped,
        genes_with_cds: cds_genes.len(),
    })
}

struct InvalidTranscript {
    gene_id: String,
    transcript_id: String,
    seqid: String,
}

fn build_transcript(
    mut builder: TranscriptBuilder,
    genes: &HashMap<String, GeneBuilder>,
) -> Result<TranscriptModel, (InvalidTranscript, String)> {
    let seqid = builder.seqid.clone().unwrap_or_default();
    let invalid = |reason: String| {
        (
            InvalidTranscript {
                gene_id: builder.gene_id.clone(),
                transcript_id: builder.transcript_id.clone(),
                seqid: seqid.clone(),
            },
            reason,
        )
    };
    let strand = match builder.strand {
        Some(value) => value,
        None => return Err(invalid("transcript has no strand".into())),
    };
    if let Some(reason) = builder.error.clone() {
        return Err(invalid(reason));
    }
    builder.cds.sort_by(|a, b| {
        if strand == Strand::Plus {
            (a.start0, a.end0).cmp(&(b.start0, b.end0))
        } else {
            (b.start0, b.end0).cmp(&(a.start0, a.end0))
        }
    });
    builder
        .cds
        .dedup_by(|a, b| a.start0 == b.start0 && a.end0 == b.end0);
    if builder.cds.iter().any(|segment| segment.seqid != seqid) {
        return Err(invalid("CDS segments use multiple sequence IDs".into()));
    }

    let cds_start = builder
        .cds
        .iter()
        .map(|segment| segment.start0)
        .min()
        .unwrap_or(0);
    let cds_end = builder
        .cds
        .iter()
        .map(|segment| segment.end0)
        .max()
        .unwrap_or(0);
    let start0 = builder.start0.unwrap_or(cds_start).min(cds_start);
    let end0 = builder.end0.unwrap_or(cds_end).max(cds_end);
    let cds_length_bp: u64 = builder.cds.iter().map(CdsSegment::len).sum();
    if cds_length_bp == 0 {
        return Err(invalid("CDS total length is zero".into()));
    }

    let first_phase = builder
        .cds
        .first()
        .and_then(|segment| segment.phase)
        .map_or(0, u64::from);
    let mut cumulative = 0u64;
    let mut phase_consistent = true;
    for (index, segment) in builder.cds.iter().enumerate() {
        if index > 0 {
            let coding_length = cumulative.saturating_sub(first_phase);
            let expected = ((3 - coding_length % 3) % 3) as u8;
            if segment.phase.is_some_and(|phase| phase != expected) {
                phase_consistent = false;
            }
        }
        cumulative += segment.len();
    }
    if first_phase > builder.cds.first().map_or(0, CdsSegment::len) {
        phase_consistent = false;
    }

    let gene = if let Some(gene) = genes.get(&builder.gene_id) {
        if gene.seqid != seqid || gene.strand != strand {
            return Err(invalid(
                "gene and transcript sequence/strand conflict".into(),
            ));
        }
        GeneAnnotation {
            gene_id: gene.gene_id.clone(),
            gene_name: gene.gene_name.clone(),
            gene_biotype: gene.gene_biotype.clone(),
            seqid: gene.seqid.clone(),
            start0: gene.start0,
            end0: gene.end0,
            strand: gene.strand,
        }
    } else {
        GeneAnnotation {
            gene_id: builder.gene_id.clone(),
            gene_name: None,
            gene_biotype: None,
            seqid: seqid.clone(),
            start0,
            end0,
            strand,
        }
    };

    Ok(TranscriptModel {
        gene,
        transcript_id: builder.transcript_id,
        transcript_name: builder.transcript_name,
        seqid,
        start0,
        end0,
        strand,
        cds: builder.cds,
        cds_length_bp,
        phase_consistent,
    })
}

fn parse_coordinate(value: &str, line_number: usize) -> io::Result<u64> {
    value
        .parse::<u64>()
        .map_err(|_| invalid_line(line_number, &format!("invalid coordinate '{value}'")))
}

fn invalid_line(line_number: usize, message: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("GFF3 line {line_number}: {message}"),
    )
}

fn first_attribute(attributes: &HashMap<String, Vec<String>>, keys: &[&str]) -> Option<String> {
    keys.iter().find_map(|key| {
        attributes
            .get(*key)
            .and_then(|values| values.first())
            .cloned()
    })
}

fn parse_attributes(text: &str) -> HashMap<String, Vec<String>> {
    let mut result = HashMap::new();
    for item in text.split(';') {
        let item = item.trim();
        if item.is_empty() {
            continue;
        }
        if let Some((key, value)) = item.split_once('=') {
            let values = if key == "Parent" {
                value
                    .split(',')
                    .map(percent_decode)
                    .filter(|value| !value.is_empty())
                    .collect()
            } else {
                vec![percent_decode(value)]
            };
            result.insert(key.to_string(), values);
        }
    }
    result
}

fn percent_decode(value: &str) -> String {
    let bytes = value.as_bytes();
    let mut output = Vec::with_capacity(bytes.len());
    let mut index = 0usize;
    while index < bytes.len() {
        if bytes[index] == b'%' && index + 2 < bytes.len() {
            if let (Some(high), Some(low)) =
                (hex_value(bytes[index + 1]), hex_value(bytes[index + 2]))
            {
                output.push((high << 4) | low);
                index += 3;
                continue;
            }
        }
        output.push(bytes[index]);
        index += 1;
    }
    String::from_utf8_lossy(&output).into_owned()
}

fn hex_value(value: u8) -> Option<u8> {
    match value {
        b'0'..=b'9' => Some(value - b'0'),
        b'a'..=b'f' => Some(value - b'a' + 10),
        b'A'..=b'F' => Some(value - b'A' + 10),
        _ => None,
    }
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|base| match base {
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
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn parse(text: &str) -> AnnotationSet {
        read_gff3_from(Cursor::new(text.as_bytes())).unwrap()
    }

    #[test]
    fn parses_attributes_and_multiple_parents() {
        let set = parse(concat!(
            "##gff-version 3\n",
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;Name=Gene%201;gene_biotype=protein_coding\n",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g1\n",
            "chr1\tsrc\tCDS\t1\t60\t.\t+\t0\tParent=t1,t2\n",
            "chr1\tsrc\tCDS\t61\t100\t.\t+\t0\tParent=t1,t2\n",
            "##FASTA\n",
            ">chr1\n",
            "AAAA\n"
        ));
        assert_eq!(set.transcripts.len(), 1);
        assert_eq!(set.genes_with_cds, 1);
        assert_eq!(set.transcripts[0].cds_length_bp, 100);
        assert_eq!(set.transcripts[0].gene.gene_name.as_deref(), Some("Gene 1"));
    }

    #[test]
    fn chooses_longest_transcript_with_stable_tie_break() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g1\n",
            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tParent=t2\n",
            "chr1\tsrc\tCDS\t1\t60\t.\t+\t0\tParent=t1\n"
        ));
        assert_eq!(set.transcripts[0].transcript_id, "t1");
    }

    #[test]
    fn extracts_and_reverse_complements_cds() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t1\t9\t.\t-\t.\tID=g1\n",
            "chr1\tsrc\tmRNA\t1\t9\t.\t-\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t7\t9\t.\t-\t0\tParent=t1\n",
            "chr1\tsrc\tCDS\t1\t3\t.\t-\t0\tParent=t1\n"
        ));
        let cds = set.transcripts[0].extract_cds("AAACCCGGG").unwrap();
        assert_eq!(cds, "CCCTTT");
    }

    #[test]
    fn records_phase_inconsistency() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=g1\n",
            "chr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t1\t4\t.\t+\t0\tParent=t1\n",
            "chr1\tsrc\tCDS\t5\t10\t.\t+\t1\tParent=t1\n"
        ));
        assert!(!set.transcripts[0].phase_consistent);
    }

    #[test]
    fn accepts_consistent_partial_cds_and_exposes_phase_offset() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1\n",
            "chr1\tsrc\tmRNA\t1\t10\t.\t+\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t1\t4\t.\t+\t1\tParent=t1\n",
            "chr1\tsrc\tCDS\t5\t10\t.\t+\t0\tParent=t1\n"
        ));
        assert!(set.transcripts[0].phase_consistent);
        assert_eq!(set.transcripts[0].translation_offset(), 1);
    }

    #[test]
    fn reverse_strand_partial_phase_uses_transcript_order() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t1\t10\t.\t-\t.\tID=g1\n",
            "chr1\tsrc\tmRNA\t1\t10\t.\t-\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t7\t10\t.\t-\t2\tParent=t1\n",
            "chr1\tsrc\tCDS\t1\t6\t.\t-\t1\tParent=t1\n"
        ));
        assert_eq!(set.transcripts[0].translation_offset(), 2);
        assert!(set.transcripts[0].phase_consistent);
    }

    #[test]
    fn merges_multipart_gene_records_with_same_id() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t100\t120\t.\t-\t.\tID=g1;Name=gene1\n",
            "chr1\tsrc\tgene\t140\t160\t.\t-\t.\tID=g1;Name=gene1\n",
            "chr1\tsrc\tmRNA\t100\t160\t.\t-\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t100\t160\t.\t-\t0\tParent=t1\n"
        ));
        assert_eq!(set.transcripts[0].gene.start0, 99);
        assert_eq!(set.transcripts[0].gene.end0, 160);
    }

    #[test]
    fn reads_gzip_gff3_by_magic_number() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let path = std::env::temp_dir().join(format!(
            "fastnlr-gff3-{}-{}.gff3.gz",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let text = concat!(
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n",
            "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
        );
        let file = std::fs::File::create(&path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(text.as_bytes()).unwrap();
        encoder.finish().unwrap();

        let set = read_gff3(&path).unwrap();
        std::fs::remove_file(path).unwrap();
        assert_eq!(set.transcripts.len(), 1);
        assert_eq!(set.transcripts[0].cds_length_bp, 9);
    }

    #[test]
    fn skips_unoriented_features() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=g1\n",
            "chr1\tsrc\tgene\t30\t40\t.\t?\t.\tID=g2\n",
            "chr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t1\t20\t.\t+\t0\tParent=t1\n"
        ));
        assert_eq!(set.transcripts.len(), 1);
        assert_eq!(set.transcripts[0].gene.gene_id, "g1");
    }

    #[test]
    fn skips_trans_spliced_transcript_with_mixed_cds_strands() {
        let set = parse(concat!(
            "chr1\tsrc\tgene\t1\t60\t.\t-\t.\tID=g1\n",
            "chr1\tsrc\tgene\t20\t40\t.\t+\t.\tID=g1\n",
            "chr1\tsrc\tmRNA\t1\t60\t.\t?\t.\tID=t1;Parent=g1\n",
            "chr1\tsrc\tCDS\t1\t10\t.\t-\t0\tParent=t1\n",
            "chr1\tsrc\tCDS\t20\t30\t.\t+\t0\tParent=t1\n",
            "chr1\tsrc\tgene\t100\t150\t.\t+\t.\tID=g2\n",
            "chr1\tsrc\tmRNA\t100\t150\t.\t+\t.\tID=t2;Parent=g2\n",
            "chr1\tsrc\tCDS\t100\t150\t.\t+\t0\tParent=t2\n"
        ));
        assert_eq!(set.transcripts.len(), 1);
        assert_eq!(set.transcripts[0].gene.gene_id, "g2");
        assert_eq!(set.skipped.len(), 1);
        assert!(set.skipped[0].reason.contains("conflicting strands"));
    }
}
