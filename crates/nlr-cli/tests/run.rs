//! nlr-cli run() orchestration integration test: full pipeline (chop -> translate -> scan -> assemble).

use std::io::Write;

/// Precise back-translation (to construct DNA that hits motifs).
fn back_translate(prot: &str) -> String {
    let mut map = std::collections::HashMap::new();
    map.insert('P', "CCG");
    map.insert('I', "ATT");
    map.insert('W', "TGG");
    map.insert('G', "GGT");
    map.insert('M', "ATG");
    map.insert('V', "GTT");
    map.insert('K', "AAA");
    map.insert('T', "ACT");
    map.insert('L', "CTT");
    map.insert('A', "GCT");
    map.insert('R', "CGT");
    map.insert('Y', "TAT");
    map.insert('N', "AAT");
    map.insert('D', "GAT");
    map.insert('H', "CAT");
    map.insert('F', "TTT");
    map.insert('C', "TGT");
    map.insert('S', "TCT");
    map.insert('Q', "CAG");
    map.insert('E', "GAA");
    prot.chars()
        .map(|a| map.get(&a).copied().unwrap_or("GGG"))
        .collect::<String>()
}

#[test]
fn run_produces_motifs_from_dna() {
    // Build a protein containing motif_1(21aa) + motif_6(29aa), back-translate to DNA,
    // add flanking to reach sufficient length.
    let m1 = "PIWGMGGVGKTTLARAVYNDP";
    let m6 = "HFDCRAWVCVSQQYDMKKVLRDIIQQVGG";
    let prot = format!("{}GGGGGGGGGG{}{}GGGGGGGGGG{}", m1, m6, m1, m6);
    let prot = prot + &"A".repeat(100);
    let dna = back_translate(&prot);

    // Write a temp FASTA (fragment_length defaults to 20000; a sequence < 20000 still yields one fragment).
    let mut fasta_path = std::env::temp_dir();
    fasta_path.push(format!("nlr_test_{}.fasta", std::process::id()));
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1").unwrap();
        for chunk in dna.as_bytes().chunks(80) {
            writeln!(f, "{}", String::from_utf8_lossy(chunk)).unwrap();
        }
    }

    // mot/store = None -> uses built-in embedded config.
    let config = nlr_cli::RunConfig::new(fasta_path.clone(), None, None);

    let result = nlr_cli::run(&config).expect("run should succeed");

    eprintln!(
        "motif sequence count: {}, NLR count: {}",
        result.motifs_by_seq.len(),
        result.nlrs.len()
    );
    for (seq, motifs) in &result.motifs_by_seq {
        let ids: Vec<u8> = motifs.iter().map(|m| m.id).collect();
        eprintln!("seq {}: {} motifs, ids={:?}", seq, motifs.len(), ids);
    }

    // Should scan at least motif_1.
    let total: usize = result.motifs_by_seq.values().map(|v| v.len()).sum();
    assert!(total > 0, "should scan at least one motif");

    let _ = std::fs::remove_file(&fasta_path);
}

/// run()/run_tracked 返回的位点必须是收尾精修后的最终结果（HMM N 端补充 + 高精度
/// 过滤 + 重命名），服务层/桌面端不经过 write_genome_outputs 也应拿到一致口径；
/// 同时验证 finalize_result 的幂等性（重复调用为零操作）。
#[test]
fn run_returns_finalized_loci_and_finalize_is_idempotent() {
    let m1 = "PIWGMGGVGKTTLARAVYNDP";
    let m6 = "HFDCRAWVCVSQQYDMKKVLRDIIQQVGG";
    let prot = format!("{}GGGGGGGGGG{}{}GGGGGGGGGG{}", m1, m6, m1, m6);
    let prot = prot + &"A".repeat(100);
    let dna = back_translate(&prot);

    let mut fasta_path = std::env::temp_dir();
    fasta_path.push(format!("nlr_finalized_{}.fasta", std::process::id()));
    {
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">chr1").unwrap();
        for chunk in dna.as_bytes().chunks(80) {
            writeln!(f, "{}", String::from_utf8_lossy(chunk)).unwrap();
        }
    }

    let config = nlr_cli::RunConfig::new(fasta_path.clone(), None, None);
    let mut result = nlr_cli::run(&config).expect("run should succeed");

    assert!(result.finalized, "run() 必须返回已完成收尾精修的位点");
    // 重命名收尾：位点名称必须与 {seq}_nlr{N} 连续编号一致。
    for (i, list) in result.nlrs.iter().enumerate() {
        assert_eq!(list.name, format!("chr1_nlr{}", i + 1));
    }

    // 幂等：对已收尾结果重复调用 finalize_result 必须是零操作。
    let nlrs_before = result.nlrs.clone();
    let annotated = nlr_cli::finalize_result(&fasta_path, &mut result, None).unwrap();
    assert_eq!(annotated, 0, "已收尾结果不应再次执行 HMM 补充");
    assert!(result.finalized);
    assert_eq!(result.nlrs.len(), nlrs_before.len(), "重复收尾不得增删位点");
    for (a, b) in result.nlrs.iter().zip(&nlrs_before) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.motifs.len(), b.motifs.len());
    }

    let _ = std::fs::remove_file(&fasta_path);
}
