//! 端到端输出回归：用真实二进制跑合成基因组，校验输出坐标与块结构自洽。
//!
//! 夹具由内置 NB-ARC 共识序列拼成，相邻 motif 的 dna_start 间距 45~87 bp，
//! 远小于播种距离阈值 500 bp，因此必然形成互相重叠的种子 —— 这正是历史上
//! 位点内 motif 重复（BED 块重复/重叠、NB-ARC 比对错位）的触发条件。

use std::process::Command;

use nlr_core::signature_def::AnnotatorSignatureDefinition;

/// 每个氨基酸固定一个密码子，覆盖全部 20 种残基；未覆盖时直接报错而不是静默兜底。
fn codon(aa: char) -> &'static str {
    match aa {
        'A' => "GCT",
        'R' => "CGT",
        'N' => "AAT",
        'D' => "GAT",
        'C' => "TGT",
        'Q' => "CAA",
        'E' => "GAA",
        'G' => "GGT",
        'H' => "CAT",
        'I' => "ATT",
        'L' => "CTT",
        'K' => "AAA",
        'M' => "ATG",
        'F' => "TTT",
        'P' => "CCG",
        'S' => "TCT",
        'T' => "ACT",
        'W' => "TGG",
        'Y' => "TAT",
        'V' => "GTT",
        other => panic!("no codon for residue {other}"),
    }
}

#[test]
fn outputs_have_no_duplicate_bed_blocks() {
    let def = AnnotatorSignatureDefinition::new();
    let protein: String = [1u8, 6, 4, 5, 10, 3, 12, 2]
        .iter()
        .map(|&id| def.consensus(id))
        .collect();
    let dna: String = protein.chars().map(codon).collect();

    let tmp = tempfile::TempDir::new().unwrap();
    let fasta = tmp.path().join("fixture.fasta");
    std::fs::write(&fasta, format!(">chr1\n{dna}\n")).unwrap();
    let out_dir = tmp.path().join("out");

    let status = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "-i",
            fasta.to_str().unwrap(),
            "-o",
            out_dir.to_str().unwrap(),
            "-t",
            "2",
        ])
        .status()
        .unwrap();
    assert!(status.success(), "fastnlr 运行失败: {status}");

    // BED：一个位点、8 个块、块之间不重叠且按坐标升序。
    let bed = std::fs::read_to_string(out_dir.join("fixture.nlr.bed")).unwrap();
    let rows: Vec<&str> = bed.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(rows.len(), 1, "夹具应只产生一个位点");

    let cols: Vec<&str> = rows[0].split('\t').collect();
    let sizes: Vec<u64> = cols[10]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect();
    let starts: Vec<u64> = cols[11]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(sizes.len(), starts.len());
    assert_eq!(cols[9].parse::<usize>().unwrap(), starts.len(), "blockCount 与块数不符");

    let mut blocks: Vec<(u64, u64)> = starts
        .iter()
        .zip(&sizes)
        .map(|(&s, &z)| (s, s + z))
        .collect();
    blocks.sort();
    for w in blocks.windows(2) {
        assert!(
            w[0].1 <= w[1].0,
            "BED 块重叠：{:?} 与 {:?}（位点内 motif 重复）",
            w[0],
            w[1]
        );
    }
    assert_eq!(blocks.len(), 8, "8 个 NB-ARC motif 应对应 8 个块");

    // txt 的 motif 列表长度必须与 BED 块数一致（两者同源）。
    let txt = std::fs::read_to_string(out_dir.join("fixture.nlr.txt")).unwrap();
    let cols: Vec<&str> = txt.trim_end().split('\t').collect();
    assert_eq!(cols[6].split(',').count(), blocks.len(), "txt motif 列表长度 != BED 块数");
}
