//! nlr-domain — 蛋白模式的结构域扫描（HMMER3）与 NLR 结构域架构判定。
//!
//! 设计要点：
//! - **结构域证据**：用 `hmmer-pure-rs`（HMMER 3.4 的纯 Rust 移植）扫描内置或用户提供的 Pfam
//!   HMM，按 Pfam GA（gathering）位分阈值判定模型是否成立，同时给出 bitscore 与 p 值。
//! - **LRR 证据**：由调用方用现有 PWM motif（motif_9/11）判定后以计数传入。实测 Pfam 的 7 个
//!   LRR 模型联合灵敏度仅 47.5%，而 PWM 在同一金标准集上为 99.5%，故 LRR 不走 HMM。
//! - **无基因组坐标**：蛋白模式没有坐标概念，输出为结构域表与架构表，不产生 GFF/BED/loci。
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

/// 内置精选模型：(accession, 名称, 类别, HMM 文本)。
///
/// 取自 InterPro/Pfam 的 HMM（见 `data/README.md` 与 `data/fetch_hmms.py`）：
/// NB-ARC 与 TIR 是 NLR 的核心结构域；Rx_N 覆盖 CNL 型 CC，RPW8 覆盖 RNL/helper 型 CC。
pub const EMBEDDED_MODELS: [(&str, &str, DomainCategory, &str); 4] = [
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
    /// N 端结构域：`CC` / `TIR` / `NB-only` / `None`。
    pub nterm: String,
    pub has_nbarc: bool,
    pub has_tir: bool,
    pub has_lrr: bool,
    /// NB-ARC + LRR 视为完整（蛋白模式没有 P-loop motif 锚点）。
    pub complete: bool,
    /// 是否算 NLR 候选：含 NB-ARC，或含 TIR 且含 LRR。
    pub is_nlr: bool,
}

/// 已知 accession 到结构域类别的映射（用户自带 HMM 用；未知返回 `NA`）。
pub fn category_for_accession(acc: &str) -> DomainCategory {
    match acc.split('.').next().unwrap_or(acc) {
        "PF00931" => DomainCategory::Nbarc,
        "PF01582" => DomainCategory::Tir,
        "PF18052" | "PF05659" => DomainCategory::Cc,
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
        let template = &self.models[0].hmm;
        let abc = Alphabet::new(template.abc_type);
        let bg = Bg::new(&abc);
        let mut pli = Pipeline::new();
        // 单序列搜索：数据库大小取 1，使 E 值等于 p 值；判定仍以 GA 位分为准。
        pli.z = 1.0;
        pli.domz = 1.0;
        pli.e_value_threshold = 10.0;
        pli.dom_e_value_threshold = 10.0;

        let prepared = self
            .models
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
/// - `lrr_hits` 是调用方用 PWM motif 得到的 LRR 命中数（>=1 即视为含 LRR）；
/// - 类别串按位置排序、折叠相邻重复，风格与 DNA 模式的 `domain_string` 一致。
pub fn classify(hmm_hits: &[DomainHit], lrr_hits: usize) -> NlrArchitecture {
    let mut passed: Vec<&DomainHit> = hmm_hits.iter().filter(|h| h.passed).collect();
    passed.sort_by(|a, b| (a.start, a.end).cmp(&(b.start, b.end)));

    let has_nbarc = passed.iter().any(|h| h.category == DomainCategory::Nbarc);
    let has_tir = passed.iter().any(|h| h.category == DomainCategory::Tir);
    let has_lrr = lrr_hits > 0;
    let first_nbarc = passed
        .iter()
        .find(|h| h.category == DomainCategory::Nbarc)
        .map(|h| h.start);

    // N 端结构域：只认位于首个 NB-ARC 之前（或没有 NB-ARC 时）的 TIR/CC，取最靠前者。
    let nterm = passed
        .iter()
        .filter(|h| matches!(h.category, DomainCategory::Tir | DomainCategory::Cc))
        .filter(|h| match first_nbarc {
            Some(nb) => h.start < nb,
            None => true,
        })
        .min_by_key(|h| h.start)
        .map(|h| h.category.as_str().to_string())
        .unwrap_or_else(|| {
            if has_nbarc {
                "NB-only".to_string()
            } else {
                "None".to_string()
            }
        });

    // 类别串：按位置顺序，LRR 作为 C 端证据追加。
    let mut parts: Vec<&'static str> = Vec::new();
    for h in &passed {
        if h.category == DomainCategory::Na {
            continue;
        }
        let s = h.category.as_str();
        if parts.last() != Some(&s) {
            parts.push(s);
        }
    }
    if has_lrr && parts.last() != Some(&"LRR") {
        parts.push("LRR");
    }

    NlrArchitecture {
        class: parts.join("-"),
        nterm,
        has_nbarc,
        has_tir,
        has_lrr,
        complete: has_nbarc && has_lrr,
        is_nlr: has_nbarc || (has_tir && has_lrr),
    }
}
