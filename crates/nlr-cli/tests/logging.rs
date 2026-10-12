//! CLI 行为回归测试：日志时间戳格式与柔性播种的默认档位。

use std::process::Command;

/// 写入一条最小 FASTA，返回 (临时目录, 输入路径, 输出路径)。
fn tiny_input() -> (tempfile::TempDir, String, String) {
    let tmp = tempfile::TempDir::new().unwrap();
    let fasta = tmp.path().join("tiny.fa");
    std::fs::write(&fasta, ">chr1\nATGGCTAGCTAGCTAGCTAGCTAGCTAGCTAA\n").unwrap();
    let out = tmp.path().join("out");
    (
        tmp,
        fasta.to_str().unwrap().to_string(),
        out.to_str().unwrap().to_string(),
    )
}

/// 运行 fastnlr 并返回 stdout（日志默认写 stdout）。
fn run(args: &[&str]) -> (bool, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
        .args(args)
        .output()
        .unwrap();
    (
        output.status.success(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
    )
}

/// 判断字符串是否为 `YYYY-MM-DD HH:MM:SS`（长度固定 19，不含小数秒）。
fn is_local_timestamp(s: &str) -> bool {
    let b = s.as_bytes();
    if b.len() != 19 {
        return false;
    }
    let d = |i: usize| b[i].is_ascii_digit();
    d(0) && d(1)
        && d(2)
        && d(3)
        && b[4] == b'-'
        && d(5)
        && d(6)
        && b[7] == b'-'
        && d(8)
        && d(9)
        && b[10] == b' '
        && d(11)
        && d(12)
        && b[13] == b':'
        && d(14)
        && d(15)
        && b[16] == b':'
        && d(17)
        && d(18)
}

#[test]
fn log_timestamp_is_local_time_without_fractional_seconds() {
    let (_tmp, fasta, out) = tiny_input();
    let (ok, log) = run(&["-i", &fasta, "-o", &out, "-t", "1", "--progress", "off"]);
    assert!(ok, "运行失败");

    let first = log
        .lines()
        .find(|l| !l.trim().is_empty())
        .expect("应有日志输出");
    let stamp = &first[..first.len().min(19)];
    assert!(
        is_local_timestamp(stamp),
        "时间戳应为 YYYY-MM-DD HH:MM:SS 且不保留小数秒，实际首行: {first}"
    );
    // 时间戳之后必须是分隔空格，确保没有 ".xxxxxx" 或 ISO 的 "T…Z" 形式。
    assert_eq!(
        first.as_bytes()[19],
        b' ',
        "时间戳与日志级别之间应为空格（不应带小数秒）: {first}"
    );
}

#[test]
fn flexible_seeding_is_enabled_by_default() {
    let (_tmp, fasta, out) = tiny_input();
    let (ok, log) = run(&["-i", &fasta, "-o", &out, "-t", "1", "--progress", "off"]);
    assert!(ok, "运行失败");
    assert!(
        log.contains("seeding: flexible (min_nbarc=3)"),
        "默认应开启柔性播种，实际日志: {log}"
    );

    // --no-flexible-seed 回退到仅精确组合。
    let (_tmp2, fasta2, out2) = tiny_input();
    let (ok, log) = run(&[
        "-i",
        &fasta2,
        "-o",
        &out2,
        "-t",
        "1",
        "--progress",
        "off",
        "--no-flexible-seed",
    ]);
    assert!(ok, "运行失败");
    assert!(
        log.contains("seeding: exact seed combinations only"),
        "--no-flexible-seed 应回退到仅精确组合，实际日志: {log}"
    );
}

#[test]
fn conflicting_seeding_flags_are_rejected() {
    let (_tmp, fasta, out) = tiny_input();
    let cases: [(&[&str], &str); 3] = [
        (
            &["--flexible-seed", "2", "--no-flexible-seed"],
            "--flexible-seed 与 --no-flexible-seed 互斥",
        ),
        (
            &["--no-flexible-seed", "--require-ploop"],
            "--require-ploop 只在柔性播种开启时有意义",
        ),
        (
            &["--relaxed-seed", "3"],
            "--relaxed-seed 与默认开启的柔性播种互斥",
        ),
    ];
    for (extra, expect) in cases {
        let mut args = vec!["-i", fasta.as_str(), "-o", out.as_str()];
        args.extend_from_slice(extra);
        let output = Command::new(env!("CARGO_BIN_EXE_fastnlr"))
            .args(&args)
            .output()
            .unwrap();
        assert!(!output.status.success(), "{extra:?} 应被拒绝");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains(expect),
            "{extra:?} 的错误信息应包含 {expect}，实际: {stderr}"
        );
    }
}
