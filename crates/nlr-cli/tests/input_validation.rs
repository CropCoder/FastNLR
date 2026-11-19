//! CLI 回归测试：非法 FASTA 必须在扫描前失败。
//!
//! Author: Jiwen Zhao (https://github.com/CropCoder)

use std::process::Command;

#[test]
fn cli_rejects_duplicate_fasta_ids_without_writing_results() {
    let tmp = tempfile::TempDir::new().unwrap();
    let fasta = tmp.path().join("duplicate.fa");
    let output_dir = tmp.path().join("out");
    std::fs::write(&fasta, ">chr1\nACGTACGT\n>chr1\nTGCATGCA\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args([
            "-i",
            fasta.to_str().unwrap(),
            "-o",
            output_dir.to_str().unwrap(),
            "-t",
            "1",
            "--progress",
            "off",
        ])
        .output()
        .unwrap();

    assert!(!output.status.success(), "duplicate IDs must fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("duplicate FASTA identifier 'chr1'"),
        "unexpected error: {stderr}"
    );
    assert!(
        !output_dir.join("duplicate.nlr.txt").exists(),
        "malformed input must not produce locus output"
    );
}
