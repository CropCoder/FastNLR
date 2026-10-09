//! nlr-domain 集成测试：用真实 NLR 蛋白（UniProt，见 `tests/data/README.md`）验证结构域扫描，
//! 并用打乱序列作为阴性对照。
//!
//! 断言里锁定了 bitscore（±0.5 bit）：`hmmer-pure-rs` 升级若改变打分，测试必须显式更新。

use nlr_domain::{classify, DomainCategory, DomainHit, HmmLibrary};

fn read_fasta(path: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq = String::new();
    for line in std::fs::read_to_string(path).unwrap().lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if !name.is_empty() {
                out.push((name.clone(), seq.clone()));
                seq.clear();
            }
            name = rest
                .split('|')
                .nth(2)
                .unwrap_or(rest)
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !name.is_empty() {
        out.push((name, seq));
    }
    out
}

fn fixtures() -> Vec<(String, String)> {
    read_fasta(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/nlr_proteins.fa"))
}

/// 确定性打乱（xorshift + Fisher-Yates），用于阴性对照。
fn shuffled(seq: &str, seed: u64) -> String {
    let mut bytes = seq.as_bytes().to_vec();
    let mut state = seed | 1;
    for i in (1..bytes.len()).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let j = (state % (i as u64 + 1)) as usize;
        bytes.swap(i, j);
    }
    String::from_utf8(bytes).unwrap()
}

#[test]
fn embedded_models_find_expected_domain_architecture() {
    let library = HmmLibrary::from_embedded();
    assert_eq!(library.models().len(), 4, "内置应为 NB-ARC/TIR/Rx_N/RPW8 四个模型");

    let mut scanner = library.scanner();
    let mut seen = 0;
    for (id, seq) in fixtures() {
        let hits = scanner.scan_protein(&id, &seq);
        let passed: Vec<_> = hits.iter().filter(|h| h.passed).collect();
        // 只算 HMM 证据时的架构（LRR 由 CLI 用 PWM 补充）。
        let arch = classify(&hits, 0);
        seen += 1;

        match id.as_str() {
            "RPM1_ARATH" => {
                let cc = passed
                    .iter()
                    .find(|h| h.model_acc == "PF18052")
                    .expect("RPM1 应有 Rx_N(CC)");
                let nb = passed
                    .iter()
                    .find(|h| h.category == DomainCategory::Nbarc)
                    .expect("RPM1 应有 NB-ARC");
                assert!((cc.bitscore - 56.2).abs() < 0.5, "Rx_N bitscore 漂移: {}", cc.bitscore);
                assert!((nb.bitscore - 193.2).abs() < 0.5, "NB-ARC bitscore 漂移: {}", nb.bitscore);
                assert_eq!(arch.nterm, "CC");
                assert_eq!(arch.class, "CC-NBARC");
                assert!(arch.has_nbarc && !arch.has_lrr && !arch.complete);
            }
            "RPS4L_ARATH" => {
                let tir = passed
                    .iter()
                    .find(|h| h.category == DomainCategory::Tir)
                    .expect("RPS4 应有 TIR");
                assert!((tir.bitscore - 166.2).abs() < 0.5, "TIR bitscore 漂移: {}", tir.bitscore);
                assert_eq!(arch.nterm, "TIR");
                assert_eq!(arch.class, "TIR-NBARC");
                assert!(arch.has_tir && arch.has_nbarc && !arch.complete);
            }
            "ADR1_ARATH" => {
                let rpw8 = passed
                    .iter()
                    .find(|h| h.model_acc == "PF05659")
                    .expect("ADR1 应有 RPW8");
                assert!((rpw8.bitscore - 41.5).abs() < 0.5, "RPW8 bitscore 漂移: {}", rpw8.bitscore);
                assert_eq!(arch.nterm, "CC");
                assert_eq!(arch.class, "CC-NBARC");
            }
            other => panic!("unexpected fixture {other}"),
        }
    }
    assert_eq!(seen, 3);
}

#[test]
fn shuffled_sequences_are_not_called() {
    let library = HmmLibrary::from_embedded();
    let mut scanner = library.scanner();
    for (id, seq) in fixtures() {
        let scrambled = shuffled(&seq, 0x9E37_79B9_7F4A_7C15);
        assert_eq!(scrambled.len(), seq.len());
        let hits = scanner.scan_protein(&id, &scrambled);
        let passed = hits.iter().filter(|h| h.passed).count();
        let arch = classify(&hits, 1);
        assert_eq!(passed, 0, "{id} 打乱序列不应有通过 GA 的结构域: {hits:?}");
        assert!(!arch.is_nlr, "{id} 打乱序列不应判为 NLR");
    }
}

#[test]
fn classification_rules() {
    let mk = |cat: DomainCategory, start: usize, end: usize, passed: bool| DomainHit {
        protein_id: "p".into(),
        model_acc: "PF00000".into(),
        model_name: "m".into(),
        category: cat,
        start,
        end,
        bitscore: 50.0,
        pvalue: 1e-20,
        passed,
    };

    // CC-NBARC-LRR
    let hits = vec![
        mk(DomainCategory::Cc, 10, 100, true),
        mk(DomainCategory::Nbarc, 150, 350, true),
    ];
    let a = classify(&hits, 3);
    assert_eq!(a.class, "CC-NBARC-LRR");
    assert_eq!(a.nterm, "CC");
    assert!(a.is_nlr && a.complete);

    // TIR-NBARC（无 LRR）
    let hits = vec![
        mk(DomainCategory::Tir, 20, 190, true),
        mk(DomainCategory::Nbarc, 210, 380, true),
    ];
    let a = classify(&hits, 0);
    assert_eq!(a.class, "TIR-NBARC");
    assert_eq!(a.nterm, "TIR");
    assert!(!a.complete && a.is_nlr);

    // NBARC-LRR（无 N 端结构域）
    let a = classify(&[mk(DomainCategory::Nbarc, 100, 300, true)], 2);
    assert_eq!(a.class, "NBARC-LRR");
    assert_eq!(a.nterm, "NB-only");
    assert!(a.complete && a.is_nlr);

    // TIR-only + LRR 也算候选（TIR-only NLR）
    let a = classify(&[mk(DomainCategory::Tir, 30, 200, true)], 1);
    assert_eq!(a.class, "TIR-LRR");
    assert_eq!(a.nterm, "TIR");
    assert!(a.is_nlr && !a.complete);

    // 无域
    let a = classify(&[], 0);
    assert_eq!(a.class, "");
    assert_eq!(a.nterm, "None");
    assert!(!a.is_nlr);

    // 未通过 GA 的命中不参与判定
    let a = classify(&[mk(DomainCategory::Nbarc, 100, 300, false)], 0);
    assert!(!a.is_nlr && a.class.is_empty());

    // CC 位于 NB-ARC 之后时不计为 N 端结构域
    let hits = vec![
        mk(DomainCategory::Nbarc, 100, 300, true),
        mk(DomainCategory::Cc, 320, 420, true),
    ];
    let a = classify(&hits, 0);
    assert_eq!(a.nterm, "NB-only");
    assert_eq!(a.class, "NBARC-CC");
}

#[test]
fn user_models_override_embedded_models_by_accession() {
    let mut library = HmmLibrary::from_embedded();
    let before = library.models().len();
    // 用内置 NB-ARC 的 HMM 文本模拟“用户提供同 accession 的模型”：应覆盖而非新增。
    let added = library
        .add_text(include_str!("../data/PF00931.hmm"), "test")
        .unwrap();
    assert_eq!(added, 1);
    assert_eq!(library.models().len(), before, "同 accession 应覆盖");
    // 分类映射：未知 accession 归为 NA，只报告不参与判定。
    assert_eq!(nlr_domain::category_for_accession("PF99999"), DomainCategory::Na);
}
