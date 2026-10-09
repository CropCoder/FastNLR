//! 蛋白模式端到端测试：用真实二进制跑 `--protein`，校验输出集合与参数互斥。

use std::process::Command;

/// 三个真实 NLR（CNL/TNL/RNL）组成的 fixture，与 nlr-domain 的集成测试共用。
fn fixture() -> &'static str {
    concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../nlr-domain/tests/data/nlr_proteins.fa"
    )
}

#[test]
fn protein_mode_writes_its_own_output_set() {
    let tmp = tempfile::TempDir::new().unwrap();
    let out_dir = tmp.path().join("out");

    let status = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "--protein",
            "-i",
            fixture(),
            "-o",
            out_dir.to_str().unwrap(),
            "-t",
            "2",
        ])
        .status()
        .unwrap();
    assert!(status.success(), "protein mode 运行失败: {status}");

    // 蛋白模式的 4 类输出。
    for name in ["nlr_proteins.domains.tsv", "nlr_proteins.nlr.tsv", "nlr_proteins.summary.txt"] {
        assert!(out_dir.join(name).is_file(), "缺少输出文件 {name}");
    }
    for name in ["01-nlr-types.png", "02-domain-counts.png"] {
        assert!(
            out_dir.join("nlr_proteins.plots").join(name).is_file(),
            "缺少统计图 {name}"
        );
    }
    // 不应产生基因组坐标类输出。
    for name in ["nlr_proteins.nlr.gff", "nlr_proteins.nlr.bed", "nlr_proteins.motifs.bed",
                 "nlr_proteins.loci.fasta", "nlr_proteins.nbarc.fasta", "nlr_proteins.stats.tsv"] {
        assert!(!out_dir.join(name).exists(), "蛋白模式不应生成 {name}");
    }

    // nlr.tsv：三条均为 NLR 候选，类别与列数正确。
    let nlr = std::fs::read_to_string(out_dir.join("nlr_proteins.nlr.tsv")).unwrap();
    let rows: Vec<Vec<&str>> = nlr
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3, "fixture 应产生 3 条 NLR 候选");
    for r in &rows {
        assert_eq!(r.len(), 10, "nlr.tsv 应为 10 列: {r:?}");
    }
    let class_of = |needle: &str| {
        rows.iter()
            .find(|r| r[0].contains(needle))
            .map(|r| r[2])
            .unwrap_or_else(|| panic!("nlr.tsv 中找不到 {needle}"))
    };
    assert_eq!(class_of("RPM1_ARATH"), "CC-NBARC-LRR");
    assert_eq!(class_of("RPS4L_ARATH"), "TIR-NBARC-LRR");
    assert_eq!(class_of("ADR1_ARATH"), "CC-NBARC-LRR");

    // domains.tsv：RPM1 的 NB-ARC 命中应通过 GA 且坐标落在蛋白范围内。
    let dom = std::fs::read_to_string(out_dir.join("nlr_proteins.domains.tsv")).unwrap();
    let header = dom.lines().next().unwrap();
    assert_eq!(header.split('\t').count(), 10);
    let nbarc = dom
        .lines()
        .find(|l| l.contains("RPM1_ARATH") && l.contains("PF00931"))
        .expect("RPM1 应有 NB-ARC 命中");
    let cols: Vec<&str> = nbarc.split('\t').collect();
    assert_eq!(cols[1], "hmm");
    assert_eq!(cols[4], "NBARC");
    assert_eq!(cols[9], "true", "RPM1 的 NB-ARC 应通过 GA");
    let (start, end): (usize, usize) = (cols[5].parse().unwrap(), cols[6].parse().unwrap());
    assert!(start < end && end <= 926);
    // LRR 由 PWM 提供。
    assert!(dom.lines().any(|l| l.contains("RPM1_ARATH") && l.contains("\tpwm\t")));
}

#[test]
fn protein_mode_rejects_genome_only_flags() {
    let tmp = tempfile::TempDir::new().unwrap();
    let out_dir = tmp.path().join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "--protein",
            "-i",
            fixture(),
            "-o",
            out_dir.to_str().unwrap(),
            "--checkpoint",
            tmp.path().join("ckpt").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success(), "冲突参数应导致失败");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--checkpoint"), "错误信息应指出冲突参数: {stderr}");
}

#[test]
fn hmm_flag_requires_protein_mode() {
    let tmp = tempfile::TempDir::new().unwrap();
    let fasta = tmp.path().join("dna.fa");
    std::fs::write(&fasta, ">chr1\nATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAA\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "-i",
            fasta.to_str().unwrap(),
            "-o",
            tmp.path().join("out").to_str().unwrap(),
            "--hmm",
            fixture(),
        ])
        .output()
        .unwrap();
    assert!(!output.status.success(), "--hmm 在 DNA 模式下应报错");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("--hmm"), "错误信息应指出 --hmm: {stderr}");
}
