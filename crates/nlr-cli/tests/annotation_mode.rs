//! GFF3 注释模式端到端测试。

use std::path::{Path, PathBuf};
use std::process::Command;

/// 读取 `nlr-domain` 的真实 NLR fixture 第一条蛋白。
fn fixture_protein() -> String {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../nlr-domain/tests/data/nlr_proteins.fa");
    let text = std::fs::read_to_string(path).unwrap();
    let mut sequence = String::new();
    let mut in_first = false;
    for line in text.lines() {
        if line.starts_with('>') {
            if in_first {
                break;
            }
            in_first = true;
            continue;
        }
        if in_first {
            sequence.push_str(line.trim());
        }
    }
    sequence
}

fn backtranslate(protein: &str) -> String {
    protein
        .chars()
        .map(|amino_acid| match amino_acid {
            'A' => "GCT",
            'C' => "TGT",
            'D' => "GAT",
            'E' => "GAA",
            'F' => "TTT",
            'G' => "GGT",
            'H' => "CAT",
            'I' => "ATT",
            'K' => "AAA",
            'L' => "CTG",
            'M' => "ATG",
            'N' => "AAT",
            'P' => "CCT",
            'Q' => "CAA",
            'R' => "CGT",
            'S' => "TCT",
            'T' => "ACT",
            'V' => "GTT",
            'W' => "TGG",
            'Y' => "TAT",
            '*' => "TAA",
            _ => "NNN",
        })
        .collect()
}

fn write_fixture(dir: &Path) -> (PathBuf, PathBuf, usize) {
    let protein = fixture_protein();
    let dna = backtranslate(&protein);
    let codon_count = protein.len();
    let first_codons = (codon_count / 2).max(1);
    let split = first_codons * 3;
    let intron = "GGGGGGGGGGGGGGGGGGGG";
    let genome = format!("{}{}{}", &dna[..split], intron, &dna[split..]);
    let genome_path = dir.join("genome.fa");
    std::fs::write(&genome_path, format!(">chr1\n{genome}\n")).unwrap();

    let second_start = split + intron.len() + 1;
    let gff = format!(
        concat!(
            "##gff-version 3\n",
            "chr1\ttest\tgene\t1\t{genome_len}\t.\t+\t.\tID=g1;Name=RPM1;gene_biotype=protein_coding\n",
            "chr1\ttest\tmRNA\t1\t{genome_len}\t.\t+\t.\tID=t1;Parent=g1;Name=RPM1.1\n",
            "chr1\ttest\tCDS\t1\t{split}\t.\t+\t0\tParent=t1\n",
            "chr1\ttest\tCDS\t{second_start}\t{genome_len}\t.\t+\t0\tParent=t1\n"
        ),
        genome_len = genome.len(),
        split = split,
        second_start = second_start
    );
    let gff_path = dir.join("annotation.gff3");
    std::fs::write(&gff_path, gff).unwrap();
    (genome_path, gff_path, protein.len())
}

fn run_annotation(
    genome: &Path,
    gff: &Path,
    output: &Path,
    threads: usize,
    batch: usize,
) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "-i",
            genome.to_str().unwrap(),
            "--gff",
            gff.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "-t",
            &threads.to_string(),
            "-n",
            &batch.to_string(),
            "--progress",
            "off",
        ])
        .output()
        .unwrap()
}

#[test]
fn gff_mode_translates_multiexon_nlr_and_writes_complete_output_set() {
    let tmp = tempfile::TempDir::new().unwrap();
    let (genome, gff, protein_length) = write_fixture(tmp.path());
    let output = tmp.path().join("out");
    let result = run_annotation(&genome, &gff, &output, 2, 1);
    assert!(
        result.status.success(),
        "GFF3 annotation failed: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    for suffix in [
        ".gff.cds.fasta",
        ".gff.proteins.fasta",
        ".gff.domains.tsv",
        ".gff.genes.tsv",
        ".gff.nlr.tsv",
        ".gff.nlr.gff3",
        ".gff.skipped.tsv",
        ".gff.summary.txt",
        ".gff.summary.json",
    ] {
        assert!(
            output.join(format!("genome{suffix}")).is_file(),
            "missing output genome{suffix}"
        );
    }
    for dna_only in [
        ".nlr.txt",
        ".nlr.gff",
        ".nlr.bed",
        ".motifs.bed",
        ".loci.fasta",
        ".nbarc.fasta",
    ] {
        assert!(
            !output.join(format!("genome{dna_only}")).exists(),
            "GFF mode unexpectedly wrote genome{dna_only}"
        );
    }

    let genes = std::fs::read_to_string(output.join("genome.gff.genes.tsv")).unwrap();
    let row: Vec<&str> = genes
        .lines()
        .find(|line| !line.starts_with('#'))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(row.len(), 26);
    assert_eq!(row[0], "g1");
    assert_eq!(row[3], "t1");
    assert_eq!(row[11], protein_length.to_string());
    assert_eq!(row[12], "CC-NBARC-LRR");
    assert_eq!(row[13], "CNL");
    assert_eq!(row[9], "2");
    assert_eq!(row[15], "true");
    assert_eq!(row[16], "false");
    assert_eq!(row[17], "true");
    assert_eq!(row[18], "true");
    assert_eq!(row[19], "true");

    let domains = std::fs::read_to_string(output.join("genome.gff.domains.tsv")).unwrap();
    assert!(domains
        .lines()
        .any(|line| line.contains("PF00931") && line.contains("\ttrue")));
    assert!(!domains.lines().any(|line| line.contains("\tpwm\tmotif_")));
    assert!(domains
        .lines()
        .find(|line| line.contains("PF00931"))
        .is_some_and(|line| line
            .split('\t')
            .nth(12)
            .is_some_and(|value| !value.is_empty())));

    let nlr_gff = std::fs::read_to_string(output.join("genome.gff.nlr.gff3")).unwrap();
    assert!(nlr_gff.contains("\tNLR_gene\t"));
    assert!(nlr_gff.contains("Family=CNL"));
    assert!(nlr_gff.contains("HasCtermLRR=true"));
    assert!(nlr_gff.contains("Complete=true"));
}

#[test]
fn gff_mode_is_deterministic_across_threads_and_batches() {
    let tmp = tempfile::TempDir::new().unwrap();
    let (genome, gff, _) = write_fixture(tmp.path());
    let single = tmp.path().join("single");
    let parallel = tmp.path().join("parallel");
    assert!(run_annotation(&genome, &gff, &single, 1, 1)
        .status
        .success());
    assert!(run_annotation(&genome, &gff, &parallel, 4, 3)
        .status
        .success());
    for suffix in [".gff.genes.tsv", ".gff.domains.tsv", ".gff.nlr.gff3"] {
        assert_eq!(
            std::fs::read(single.join(format!("genome{suffix}"))).unwrap(),
            std::fs::read(parallel.join(format!("genome{suffix}"))).unwrap(),
            "different output for {suffix}"
        );
    }
}

#[test]
fn gff_mode_rejects_dna_only_flags_and_requires_gff_for_hmm() {
    let tmp = tempfile::TempDir::new().unwrap();
    let (genome, gff, _) = write_fixture(tmp.path());
    let checkpoint = tmp.path().join("ckpt");
    let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "-i",
            genome.to_str().unwrap(),
            "--gff",
            gff.to_str().unwrap(),
            "-o",
            tmp.path().join("out").to_str().unwrap(),
            "--checkpoint",
            checkpoint.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());

    for flag in [
        vec!["--motif-accept-p", "0.000001"],
        vec!["--motif-prelim-p", "0.01"],
    ] {
        let rejected_output = tmp.path().join("pwm-rejected");
        let mut args = vec![
            "-i",
            genome.to_str().unwrap(),
            "--gff",
            gff.to_str().unwrap(),
            "-o",
            rejected_output.to_str().unwrap(),
        ];
        args.extend(flag);
        let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("PWM"));
    }
    assert!(String::from_utf8_lossy(&output.stderr).contains("--checkpoint"));

    let hmm_without_gff = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "-i",
            genome.to_str().unwrap(),
            "-o",
            tmp.path().join("out2").to_str().unwrap(),
            "--hmm",
            gff.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!hmm_without_gff.status.success());
    assert!(String::from_utf8_lossy(&hmm_without_gff.stderr).contains("--gff"));
}

#[test]
fn direct_protein_mode_is_removed() {
    let tmp = tempfile::TempDir::new().unwrap();
    let protein = tmp.path().join("proteins.faa");
    std::fs::write(&protein, ">p1\nMKT\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "--protein",
            "-i",
            protein.to_str().unwrap(),
            "-o",
            tmp.path().join("out").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--protein"));

    let protein_without_flag = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "-i",
            protein.to_str().unwrap(),
            "-o",
            tmp.path().join("out2").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!protein_without_flag.status.success());
    assert!(String::from_utf8_lossy(&protein_without_flag.stderr).contains("genome FASTA"));
}
