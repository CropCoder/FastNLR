//! nlr-domain — 蛋白模式的结构域扫描（HMMER3）与 NLR 结构域架构判定。
//!
//! 设计要点：
//! - **结构域证据**：用 `hmmer-pure-rs`（HMMER 3.4 的纯 Rust 移植）扫描内置或用户提供的 Pfam
//!   HMM，按 Pfam GA（gathering）位分阈值判定模型是否成立，同时给出 bitscore 与 p 值。
//! - **有序结构域证据**：NB-ARC、TIR、CC 与 LRR 均来自通过 Pfam GA 的 HMM 命中。
//!   LRR 只有位于 NB-ARC（无 NB-ARC 时为 TIR）之后才作为 C 端证据参与完整 NLR 判定。
//! - **序列层坐标**：本 crate 输出蛋白上的 1-based 结构域坐标；GFF3 模式由 `nlr-cli`
//!   将坐标映射到 CDS 和基因组区间。
//!
//! `hmmer-pure-rs` 精确锁定 0.7.4（`=` 版本约束）：其得分直接参与回归测试，禁止静默升级。

use std::io::{self, BufReader, Cursor};
use std::path::Path;

use hmmer_pure_rs::hmmfile;
use hmmer_pure_rs::profile::{profile_config, P7_LOCAL};
use hmmer_pure_rs::sequence::Sequence;
use hmmer_pure_rs::{Alphabet, Bg, Hmm, OProfile, Pipeline, Profile, TopHits};

/// 复用 `nlr-core` 的结构域类别枚举，保证与 DNA 模式的类别串一致。
pub use nlr_core::signature_def::DomainCategory;

/// 内置模型库版本。
pub const MODEL_LIBRARY_VERSION: &str = "FastNLR-Pfam-HMM-23-2026-11-11";

/// 内置模型：(accession, 名称, 类别, HMM 文本)。
///
/// 取自 InterPro/Pfam 的 HMM（见 `data/README.md` 与 `data/fetch_hmms.py`）：
/// NB-ARC 与 TIR 是 NLR 的核心结构域；Rx_N 覆盖 CNL 型 CC，RPW8 覆盖 RNL/helper 型 CC。
/// LRR 重复和 LRR 端域模型只提供位置明确的 LRR 证据，分类仍使用统一的有序域规则。
pub const EMBEDDED_MODELS: [(&str, &str, DomainCategory, &str); 23] = [
    (
        "PF00931",
        "NB-ARC",
        DomainCategory::Nbarc,
        include_str!("../data/PF00931.hmm"),
    ),
    (
        "PF01582",
        "TIR",
        DomainCategory::Tir,
        include_str!("../data/PF01582.hmm"),
    ),
    (
        "PF18052",
        "Rx_N",
        DomainCategory::Cc,
        include_str!("../data/PF18052.hmm"),
    ),
    (
        "PF05659",
        "RPW8",
        DomainCategory::Cc,
        include_str!("../data/PF05659.hmm"),
    ),
    (
        "PF12061",
        "R1-type",
        DomainCategory::Cc,
        include_str!("../data/PF12061.hmm"),
    ),
    (
        "PF00560",
        "LRR_1",
        DomainCategory::Lrr,
        include_str!("../data/PF00560.hmm"),
    ),
    (
        "PF07723",
        "LRR_2",
        DomainCategory::Lrr,
        include_str!("../data/PF07723.hmm"),
    ),
    (
        "PF07725",
        "LRR_3",
        DomainCategory::Lrr,
        include_str!("../data/PF07725.hmm"),
    ),
    (
        "PF12799",
        "LRR_4",
        DomainCategory::Lrr,
        include_str!("../data/PF12799.hmm"),
    ),
    (
        "PF13516",
        "LRR_6",
        DomainCategory::Lrr,
        include_str!("../data/PF13516.hmm"),
    ),
    (
        "PF13855",
        "LRR_8",
        DomainCategory::Lrr,
        include_str!("../data/PF13855.hmm"),
    ),
    (
        "PF14580",
        "LRR_9",
        DomainCategory::Lrr,
        include_str!("../data/PF14580.hmm"),
    ),
    (
        "PF18805",
        "LRR_10",
        DomainCategory::Lrr,
        include_str!("../data/PF18805.hmm"),
    ),
    (
        "PF18831",
        "LRR_11",
        DomainCategory::Lrr,
        include_str!("../data/PF18831.hmm"),
    ),
    (
        "PF18837",
        "LRR_12",
        DomainCategory::Lrr,
        include_str!("../data/PF18837.hmm"),
    ),
    (
        "PF23598",
        "LRR_14",
        DomainCategory::Lrr,
        include_str!("../data/PF23598.hmm"),
    ),
    (
        "PF24969",
        "LRR_15",
        DomainCategory::Lrr,
        include_str!("../data/PF24969.hmm"),
    ),
    (
        "PF23247",
        "LRR_RPS2",
        DomainCategory::Lrr,
        include_str!("../data/PF23247.hmm"),
    ),
    (
        "PF23286",
        "LRR_13",
        DomainCategory::Lrr,
        include_str!("../data/PF23286.hmm"),
    ),
    (
        "PF01462",
        "LRRNT",
        DomainCategory::Lrr,
        include_str!("../data/PF01462.hmm"),
    ),
    (
        "PF08263",
        "LRRNT_2",
        DomainCategory::Lrr,
        include_str!("../data/PF08263.hmm"),
    ),
    (
        "PF01463",
        "LRRCT",
        DomainCategory::Lrr,
        include_str!("../data/PF01463.hmm"),
    ),
    (
        "PF16920",
        "LRRCT_2",
        DomainCategory::Lrr,
        include_str!("../data/PF16920.hmm"),
    ),
];

/// Pfam GA 未设置时 HMMER 写入的哨兵值（`-999999.99`）。
const CUTOFF_UNSET: f32 = -1.0e5;

/// 一个可扫描的结构域模型。
#[derive(Debug, Clone)]
pub struct DomainModel {
    pub acc: String,
    pub name: String,
    pub category: DomainCategory,
    hmm: Hmm,
}

impl DomainModel {
    /// Pfam GA 阈值 `(sequence, domain)`；domain 阈值未设置时回退到 sequence 阈值。
    ///
    /// `hmmer-pure-rs` 的 `cutoff` 数组布局为
    /// `[GA_seq, GA_dom, TC_seq, TC_dom, NC_seq, NC_dom]`（已按实际值逐项核对）。
    pub fn ga(&self) -> (f32, f32) {
        let seq = self.hmm.cutoff[0];
        let dom = if self.hmm.cutoff[1] <= CUTOFF_UNSET {
            seq
        } else {
            self.hmm.cutoff[1]
        };
        (seq, dom)
    }
}

/// 单条结构域命中（坐标为蛋白序列上的 1-based 闭区间，对应 HMMER 的 `iali..jali`）。
#[derive(Debug, Clone, PartialEq)]
pub struct DomainHit {
    pub protein_id: String,
    pub model_acc: String,
    pub model_name: String,
    pub category: DomainCategory,
    pub start: usize,
    pub end: usize,
    pub bitscore: f32,
    /// p 值 = `exp(hit.lnp)`（HMMER 的对数 p 值）。
    pub pvalue: f64,
    /// 是否通过该模型的 Pfam GA 阈值（序列级与结构域级同时满足）。
    pub passed: bool,
}

/// NLR 结构域架构判定结果。
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct NlrArchitecture {
    /// 按位置折叠相邻重复后的结构域串，如 `CC-NBARC-LRR`。
    pub class: String,
    /// 生物学家族：TNL / RNL / CNL / NB-only / TIR-only / Other。
    /// RNL 与 CNL 可共享 `CC-NBARC-LRR` 结构，但 N 端 CC 分别来自 RPW8 与 Rx_N/R1。
    pub family: String,
    /// N 端结构域：`CC` / `TIR` / `NB-only` / `None`。
    pub nterm: String,
    pub has_nbarc: bool,
    pub has_tir: bool,
    /// 任意位置存在通过 GA 的 LRR 命中（包含错误的 N 端位置）。
    pub has_lrr_any: bool,
    /// 存在位于 NB-ARC（无 NB-ARC 时位于 TIR）结构域块之后的 LRR。
    pub has_cterm_lrr: bool,
    /// 存在 NB-ARC 且存在 C 端 LRR；不等同于生物学 family。
    pub complete: bool,
    /// 是否算 NLR 候选：含 NB-ARC，或含 TIR 且含 C 端 LRR。
    pub is_nlr: bool,
}

/// 已知 accession 到结构域类别的映射（用户自带 HMM 用；未知返回 `NA`）。
pub fn category_for_accession(acc: &str) -> DomainCategory {
    match acc.split('.').next().unwrap_or(acc) {
        "PF00931" => DomainCategory::Nbarc,
        "PF01582" => DomainCategory::Tir,
        "PF18052" | "PF05659" | "PF12061" => DomainCategory::Cc,
        "PF00560" | "PF07723" | "PF07725" | "PF12799" | "PF13516" | "PF13855" | "PF14580"
        | "PF18805" | "PF18831" | "PF18837" | "PF23598" | "PF24969" | "PF23247" | "PF23286"
        | "PF01462" | "PF08263" | "PF01463" | "PF16920" => DomainCategory::Lrr,
        _ => DomainCategory::Na,
    }
}

/// 结构域模型库：内置模型 + 用户通过 `--hmm` 追加的模型。
#[derive(Debug, Clone)]
pub struct HmmLibrary {
    models: Vec<DomainModel>,
}

impl Default for HmmLibrary {
    fn default() -> Self {
        Self::from_embedded()
    }
}

impl HmmLibrary {
    /// 仅内置的精选模型（零配置路径）。
    pub fn from_embedded() -> Self {
        let models = EMBEDDED_MODELS
            .iter()
            .map(|(acc, name, category, text)| DomainModel {
                acc: (*acc).to_string(),
                name: (*name).to_string(),
                category: *category,
                hmm: parse_single_hmm(text)
                    .unwrap_or_else(|e| panic!("embedded HMM {acc} is invalid: {e}")),
            })
            .collect();
        HmmLibrary { models }
    }

    /// 追加用户 HMM 文件（可含多个模型）：同 accession 覆盖内置模型，新 accession 追加。
    /// 返回实际追加/覆盖的模型数量。
    pub fn add_file(&mut self, path: &Path) -> io::Result<usize> {
        let text = std::fs::read_to_string(path)?;
        self.add_text(&text, &path.display().to_string())
    }

    /// 从 HMM 文本追加模型（与 `add_file` 同一逻辑，便于测试与内嵌数据）。
    pub fn add_text(&mut self, text: &str, label: &str) -> io::Result<usize> {
        let hmms = parse_hmms(text)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("{label}: {e}")))?;
        let mut added = 0usize;
        for hmm in hmms {
            // 统一去掉版本后缀（PF00931.29 -> PF00931），保证与内置模型/类别映射同名同键。
            let raw = hmm.acc.clone().unwrap_or_else(|| hmm.name.clone());
            let acc = raw.split('.').next().unwrap_or(&raw).to_string();
            let model = DomainModel {
                acc: acc.clone(),
                name: hmm.name.clone(),
                category: category_for_accession(&acc),
                hmm,
            };
            match self.models.iter_mut().find(|m| m.acc == acc) {
                Some(slot) => *slot = model,
                None => self.models.push(model),
            }
            added += 1;
        }
        Ok(added)
    }

    pub fn models(&self) -> &[DomainModel] {
        &self.models
    }

    /// 构造扫描器；**必须在实际使用它的线程内调用**（`OProfile` 既不是 `Clone` 也不是 `Send`）。
    pub fn scanner(&self) -> Scanner {
        self.scanner_for(|_| true)
    }

    /// 只扫描核心结构域，避免 DNA 模式的六框 N 端补全执行高成本 LRR 模型。
    pub fn core_scanner(&self) -> Scanner {
        self.scanner_for(|model| model.category != DomainCategory::Lrr)
    }

    /// 构造两阶段蛋白扫描器：所有蛋白扫描核心域，仅 NLR 候选追加扫描 LRR。
    pub fn protein_scanner(&self) -> ProteinScanner {
        ProteinScanner {
            core: self.core_scanner(),
            lrr: self.scanner_for(|model| model.category == DomainCategory::Lrr),
        }
    }

    fn scanner_for(&self, include: impl Fn(&DomainModel) -> bool) -> Scanner {
        let selected: Vec<&DomainModel> =
            self.models.iter().filter(|model| include(model)).collect();
        let template = &selected[0].hmm;
        let abc = Alphabet::new(template.abc_type);
        let bg = Bg::new(&abc);
        let mut pli = Pipeline::new();
        // 单序列搜索：数据库大小取 1，使 E 值等于 p 值；判定仍以 GA 位分为准。
        pli.z = 1.0;
        pli.domz = 1.0;
        pli.e_value_threshold = 10.0;
        pli.dom_e_value_threshold = 10.0;

        let prepared = selected
            .iter()
            .map(|m| {
                let mut gm = Profile::new(m.hmm.m, &abc);
                profile_config(&m.hmm, &bg, &mut gm, 400, P7_LOCAL);
                let om = OProfile::convert(&gm);
                let (ga_seq, ga_dom) = m.ga();
                PreparedModel {
                    acc: m.acc.clone(),
                    name: m.name.clone(),
                    category: m.category,
                    hmm: m.hmm.clone(),
                    gm,
                    om,
                    ga_seq,
                    ga_dom,
                }
            })
            .collect();
        Scanner {
            abc,
            bg,
            pli,
            prepared,
        }
    }
}

/// 线程本地的两阶段蛋白扫描器。
pub struct ProteinScanner {
    core: Scanner,
    lrr: Scanner,
}

impl ProteinScanner {
    /// 先扫描核心域；只有含通过 GA 的 NB-ARC 或 TIR 时才扫描 LRR 模型。
    pub fn scan_protein(&mut self, protein_id: &str, seq: &str) -> Vec<DomainHit> {
        let mut hits = self.core.scan_protein(protein_id, seq);
        if hits.iter().any(|hit| {
            hit.passed && matches!(hit.category, DomainCategory::Nbarc | DomainCategory::Tir)
        }) {
            hits.extend(self.lrr.scan_protein(protein_id, seq));
            hits.sort_by(|a, b| {
                (a.start, a.end, &a.model_acc).cmp(&(b.start, b.end, &b.model_acc))
            });
        }
        hits
    }
}

/// 一个模型在线程内的运行时状态。
struct PreparedModel {
    acc: String,
    name: String,
    category: DomainCategory,
    hmm: Hmm,
    gm: Profile,
    om: OProfile,
    ga_seq: f32,
    ga_dom: f32,
}

/// 线程本地的结构域扫描器。
pub struct Scanner {
    abc: Alphabet,
    bg: Bg,
    pli: Pipeline,
    prepared: Vec<PreparedModel>,
}

impl Scanner {
    /// 扫描一条蛋白序列，返回所有被流水线报告的结构域命中（含未通过 GA 的，便于诊断）。
    pub fn scan_protein(&mut self, protein_id: &str, seq: &str) -> Vec<DomainHit> {
        let mut out = Vec::new();
        if seq.is_empty() {
            return out;
        }
        let dsq = self.abc.digitize(seq.as_bytes());
        let sq = Sequence {
            name: protein_id.into(),
            acc: String::new(),
            desc: String::new(),
            dsq,
            n: seq.len(),
            l: seq.len(),
            taxid: -1,
        };
        for m in self.prepared.iter_mut() {
            let mut th = TopHits::new();
            self.pli.new_model(&m.gm);
            self.pli
                .run(&mut m.gm, &mut m.om, &self.bg, &m.hmm, &sq, &mut th);
            for hit in &th.hits {
                let pvalue = hit.lnp.exp();
                for d in &hit.dcl {
                    if d.iali <= 0 || d.jali < d.iali {
                        continue;
                    }
                    out.push(DomainHit {
                        protein_id: protein_id.to_string(),
                        model_acc: m.acc.clone(),
                        model_name: m.name.clone(),
                        category: m.category,
                        start: d.iali as usize,
                        end: d.jali as usize,
                        bitscore: d.bitscore,
                        pvalue,
                        passed: hit.score >= m.ga_seq && d.bitscore >= m.ga_dom,
                    });
                }
            }
        }
        out.sort_by(|a, b| (a.start, a.end, &a.model_acc).cmp(&(b.start, b.end, &b.model_acc)));
        out
    }
}

/// 解析 HMM 文本中的第一个模型。
fn parse_single_hmm(text: &str) -> Result<Hmm, String> {
    let mut hmms = parse_hmms(text)?;
    if hmms.is_empty() {
        return Err("no HMM record found".to_string());
    }
    Ok(hmms.remove(0))
}

/// 解析 HMM 文本中的全部模型（自动识别 ASCII / 二进制格式）。
fn parse_hmms(text: &str) -> Result<Vec<Hmm>, String> {
    let reader = BufReader::new(Cursor::new(text.as_bytes()));
    hmmfile::read_hmms_auto(reader).map_err(|e| e.to_string())
}

/// 依据结构域证据判定 NLR 架构。
///
/// - 只有 `passed` 的 HMM 命中参与判定；
/// - 类别串与 DNA 模式共用 `nlr_core::signature_def::canonical_domain_string`（有序域模型），
///   保证同一基因在两种模式下判定口径一致。
pub fn classify(hmm_hits: &[DomainHit]) -> NlrArchitecture {
    let mut passed: Vec<&DomainHit> = hmm_hits.iter().filter(|h| h.passed).collect();
    passed.sort_by(|a, b| (a.start, a.end, &a.model_acc).cmp(&(b.start, b.end, &b.model_acc)));

    let has_nbarc = passed.iter().any(|h| h.category == DomainCategory::Nbarc);
    let has_tir = passed.iter().any(|h| h.category == DomainCategory::Tir);
    let has_lrr_any = passed.iter().any(|h| h.category == DomainCategory::Lrr);
    let has_cterm_lrr = passed_lrr_is_cterm(&passed);
    let selected_nterm_rpw8 = selected_nterm_uses_rpw8(&passed);

    let cats: Vec<DomainCategory> = passed.iter().map(|h| h.category).collect();
    let scores: Vec<f64> = passed.iter().map(|h| h.bitscore as f64).collect();
    let (nterm_cat, class) = nlr_core::signature_def::canonical_domain_string_with_scores(
        &cats,
        Some(&scores),
        nlr_core::signature_def::PROTEIN_REPEAT_MIN_BLOCK,
    );
    // N 端结构域对外仍是字符串：CC / TIR / NB-only / None。
    let nterm = match nterm_cat {
        Some(c) => c.as_str().to_string(),
        None if has_nbarc => "NB-only".to_string(),
        None => "None".to_string(),
    };
    let family = if has_nbarc && has_tir {
        "TNL"
    } else if has_nbarc && nterm == "CC" && selected_nterm_rpw8 {
        "RNL"
    } else if has_nbarc && nterm == "CC" {
        "CNL"
    } else if has_nbarc {
        "NB-only"
    } else if has_tir {
        "TIR-only"
    } else {
        "Other"
    };

    NlrArchitecture {
        class,
        family: family.to_string(),
        nterm,
        has_nbarc,
        has_tir,
        has_lrr_any,
        has_cterm_lrr,
        complete: has_nbarc && has_cterm_lrr,
        is_nlr: has_nbarc || (has_tir && has_cterm_lrr),
    }
}

/// C 端 LRR 必须位于首个 NB-ARC 块之后；若没有 NB-ARC，则以首个 TIR 块为锚点。
fn passed_lrr_is_cterm(passed: &[&DomainHit]) -> bool {
    let (anchor_category, anchor) = if let Some(index) = passed
        .iter()
        .position(|hit| hit.category == DomainCategory::Nbarc)
    {
        (DomainCategory::Nbarc, index)
    } else if let Some(index) = passed
        .iter()
        .position(|hit| hit.category == DomainCategory::Tir)
    {
        (DomainCategory::Tir, index)
    } else {
        return false;
    };

    let anchor_end = passed[anchor..]
        .iter()
        .take_while(|hit| hit.category == anchor_category)
        .map(|hit| hit.end)
        .max()
        .unwrap_or(passed[anchor].end);
    passed
        .iter()
        .skip(anchor + 1)
        .any(|hit| hit.category == DomainCategory::Lrr && hit.start > anchor_end)
}

/// Whether the N-terminal CC evidence selected by the ordered-domain model contains RPW8.
///
/// A passed RPW8 hit elsewhere (for example after NB-ARC) must not turn an Rx_N/R1-type
/// CNL into an RNL.
fn selected_nterm_uses_rpw8(passed: &[&DomainHit]) -> bool {
    let Some(first_nbarc) = passed
        .iter()
        .position(|h| h.category == DomainCategory::Nbarc)
    else {
        return false;
    };

    let mut blocks: Vec<(DomainCategory, usize, usize)> = Vec::new();
    for (idx, hit) in passed[..first_nbarc].iter().enumerate() {
        let cat = hit.category;
        if !matches!(cat, DomainCategory::Cc | DomainCategory::Tir) {
            continue;
        }
        match blocks.last_mut() {
            Some(last) if last.0 == cat => last.2 = idx + 1,
            _ => blocks.push((cat, idx, idx + 1)),
        }
    }
    if blocks.is_empty() {
        return false;
    }

    let mixed = blocks.iter().any(|b| b.0 == DomainCategory::Cc)
        && blocks.iter().any(|b| b.0 == DomainCategory::Tir);
    let selected = if mixed {
        blocks
            .iter()
            .enumerate()
            .max_by_key(|(idx, (_, start, end))| (end - start, *idx))
            .map(|(_, block)| block)
    } else {
        blocks.last()
    };
    let Some((category, start, end)) = selected else {
        return false;
    };
    *category == DomainCategory::Cc
        && passed[*start..*end]
            .iter()
            .any(|h| h.model_acc.split('.').next() == Some("PF05659"))
}
