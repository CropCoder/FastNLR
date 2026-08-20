//! nlr-cli run() orchestration integration test: full pipeline (chop -> translate -> scan -> assemble).

use std::io::Write;

/// Precise back-translation (to construct DNA that hits motifs).
fn back_translate(prot: &str) -> String {
    let mut map = std::collections::HashMap::new();
    map.insert('P', "CCG"); map.insert('I', "ATT"); map.insert('W', "TGG");
    map.insert('G', "GGT"); map.insert('M', "ATG"); map.insert('V', "GTT");
    map.insert('K', "AAA"); map.insert('T', "ACT"); map.insert('L', "CTT");
    map.insert('A', "GCT"); map.insert('R', "CGT"); map.insert('Y', "TAT");
    map.insert('N', "AAT"); map.insert('D', "GAT"); map.insert('H', "CAT");
    map.insert('F', "TTT"); map.insert('C', "TGT"); map.insert('S', "TCT");
    map.insert('Q', "CAG"); map.insert('E', "GAA");
    prot.chars().map(|a| map.get(&a).copied().unwrap_or("GGG")).collect::<String>()
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
