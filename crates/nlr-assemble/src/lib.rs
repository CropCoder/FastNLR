//! nlr-assemble — three-step NLR locus assembly: findSeeds -> mergeSeeds -> elongate.
//!
//! Three-step assembly algorithm:
//! - findSeeds: build seed combinations from adjacent motifs within
//!   `distanceWithinMotifCombination`;
//! - mergeSeeds: same strand + 4 dna_start combination differences within threshold + rank
//!   monotonic;
//! - elongate: extend at both ends by rank rules, with usedMotifs dedup.

use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::{AnnotatorSignatureDefinition, DomainCategory};
use nlr_core::strand::Strand;
use std::collections::{HashMap, HashSet};

/// Three-step assembly parameters.
#[derive(Debug, Clone, Copy)]
pub struct AssembleParams {
    /// Maximum distance between adjacent motifs within a seed (default 500).
    pub distance_within_motif_combination: u64,
    /// Elongation distance (default 2500).
    pub distance_for_elongating: u64,
    /// Seed merge distance (default 10000).
    pub distance_between_motif_combinations: u64,
    /// 宽松播种（新增）：连续命中串里 NB-ARC 类 motif ≥ 该数即播种；None = 沿用原版精确组合。
    pub relaxed_seed_min_nbarc: Option<usize>,
}

impl Default for AssembleParams {
    fn default() -> Self {
        AssembleParams {
            distance_within_motif_combination: 500,
            distance_for_elongating: 2500,
            distance_between_motif_combinations: 10000,
            relaxed_seed_min_nbarc: None,
        }
    }
}

/// Assemble NLR loci from all motifs on a chromosome.
///
/// Returns the named locus list (`{seq_id}_nlr{N}`).
pub fn assemble(
    seq_id: &str,
    mut motifs: Vec<Motif>,
    params: &AssembleParams,
    def: &AnnotatorSignatureDefinition,
) -> Vec<MotifList> {
    // Sort (forward strand ascending / reverse strand descending).
    motifs.sort();

    let seeds = find_seeds(
        &motifs,
        params.distance_within_motif_combination,
        def,
        params.relaxed_seed_min_nbarc,
    );
    let merged = merge_seeds(seeds, params.distance_between_motif_combinations, def);
    let mut pre_nlrs = elongate(merged, &motifs, params.distance_for_elongating, def);
    // 组装后处理：吸收被过度切分的冗余位点（与结构域判定无关，纯结构层面的收敛）。
    pre_nlrs = absorb_redundant_loci(pre_nlrs, params.distance_for_elongating, def);
    // RNL 精度保护：丢弃仅有“单个 CC + P-loop”的两 motif 短片段。
    // 完整/可延伸的 helper NLR 会先经过 absorb，再保留 LRR 或第三个结构 motif。
    pre_nlrs.retain(|list| !is_unresolved_cc_ploop_fragment(list, def));

    // Sort candidate NLRs before naming (forward strand ascending / reverse strand
    // descending), which determines the nlrN naming order.
    pre_nlrs.sort();

    let mut nlrs: Vec<MotifList> = Vec::new();
    let mut suffix = 1usize;
    for mut list in pre_nlrs {
        list.sort();
        // 去重：merge_seeds 合并重叠种子时会重复并入同一 motif（同一 id + 同一坐标），
        // 这里作为组装结果的唯一出口统一清理，保证位点内记录与坐标一一对应。
        list.remove_redundant_motifs();
        list.name = format!("{}_nlr{}", seq_id, suffix);
        suffix += 1;
        nlrs.push(list);
    }
    nlrs
}

/// 是否为未延伸的 `CC→P-loop` 两 motif 短片段。
///
/// 新增的 RNL/helper CC seed 可以恢复 ADR1 类基因，但单独一条 CC motif 与
/// P-loop 相邻时也容易产生 200-500 bp 的短片段。若该候选经合并/延伸后仍只有
/// 两个 motif 且没有 LRR，则属于未形成可解释 NLR 结构的残余 seed。
fn is_unresolved_cc_ploop_fragment(list: &MotifList, def: &AnnotatorSignatureDefinition) -> bool {
    if list.motifs.len() > 2 {
        return false;
    }
    let domain = list.domain_string(def);
    domain.starts_with("CC-NBARC") && !domain.contains("-LRR")
}

/// (1) Find seeds from adjacent motifs.
fn find_seeds(
    motifs: &[Motif],
    distance: u64,
    def: &AnnotatorSignatureDefinition,
    relaxed_seed_min_nbarc: Option<usize>,
) -> Vec<MotifList> {
    let mut seeds = Vec::new();
    let mut count = 0usize;

    let mut i = 0;
    while i < motifs.len() {
        let mut potential = vec![motifs[i].clone()];
        let mut s: Vec<u8> = vec![motifs[i].id];
        let is_seed_ids = |ids: &[u8]| {
            if let Some(cfg) = def.flexible_seed() {
                def.is_seed_flexible(ids, &cfg) || def.is_seed(ids)
            } else if let Some(n) = relaxed_seed_min_nbarc {
                def.is_seed_relaxed(ids, n)
            } else {
                def.is_seed(ids)
            }
        };

        // The original loop only checked seeding after appending a second motif.
        // With --flexible-seed 1, an isolated P-loop could never seed if the next
        // motif was more than the intra-seed distance away.
        if is_seed_ids(&s) {
            count += 1;
            let name = format!(
                "{}_nlr{}",
                motifs[i].dna_sequence_id.as_deref().unwrap_or(""),
                count
            );
            seeds.push(MotifList::new(name, potential));
            i += 1;
            continue;
        }

        let mut j = i + 1;
        let mut cur = &motifs[i];
        while j < motifs.len() {
            let next = &motifs[j];
            if cur.dna_start.abs_diff(next.dna_start) <= distance {
                s.push(next.id);
                potential.push(next.clone());
                cur = next;
                // 柔性播种启用时同时接受精确种子组合（并集语义）。
                //
                // 签名预过滤 `has_signature_flexible` 本来就是「精确签名 OR 柔性条件」的并集，
                // 播种必须与之一致：否则会出现「motif 簇因为命中精确签名而被保留下来，但播种
                // 阶段只认柔性条件、结果什么都组装不出来」的矛盾状态，进而丢掉默认配置能检出
                // 的位点。RefPlantNLR 的 Pb1（命中精确组合 [12,2,8]，但 NB-ARC 类 motif 只有
                // 2 个、达不到柔性门槛）就是这种情形，开启 --flexible-seed 后反而漏检。
                // 并集只会新增种子，不会移除任何原有种子，因此默认路径的行为完全不变。
                let seeded = is_seed_ids(&s);
                if seeded {
                    count += 1;
                    let name = format!(
                        "{}_nlr{}",
                        motifs[i].dna_sequence_id.as_deref().unwrap_or(""),
                        count
                    );
                    seeds.push(MotifList::new(name, potential));
                    break;
                }
                j += 1;
            } else {
                break;
            }
        }
        i += 1;
    }

    seeds.sort();
    seeds
}

/// (2) Merge compatible seeds.
fn merge_seeds(
    mut seeds: Vec<MotifList>,
    distance: u64,
    def: &AnnotatorSignatureDefinition,
) -> Vec<MotifList> {
    let mut merged: Vec<MotifList> = Vec::new();

    while !seeds.is_empty() {
        let mut current = seeds.remove(0);
        while !seeds.is_empty() {
            let next = &seeds[0];
            if current.can_be_merged_with(next, distance, def) {
                let other = seeds.remove(0);
                current.add_motifs(&other.motifs);
            } else {
                break;
            }
        }
        merged.push(current);
    }

    merged
}

/// (3) Elongate merged seeds at both ends.
///
/// The first/last motif is read dynamically while adding motifs, because each addition
/// re-sorts the list and shifts the elongation baseline.
///
/// 吸收被过度切分的冗余位点（组装层后处理，与结构域判定解耦）。
///
/// `merge_seeds`（rank 单调约束）与 `elongate`（只允许 rank 递增）都无法把**串联重复结构**
/// 的后半段并回主位点：例如 RPP2A 的 `TIR-NBARC-TIR-NBARC-LRR`，第二个 TIR 的 rank（1~3）
/// 低于前一段 NB-ARC 的 rank（4~12），合并必然被拒；Prv 的 `TIR-NBARC-LRR-NBARC` 同理
/// （LRR rank 14 → NB-ARC rank 4）。残留的 motif 会各自组装成小位点，表现为过度切分。
///
/// 这里按结构关系补回两类冗余位点：
/// 1. **完全嵌套**：位点 B 的跨度完全落在同链位点 A 之内。同链嵌套的两个独立 NLR 在自然界
///    不成立，这一定是同一基因的重复单元未合并（RPP2A 即此形态）。
/// 2. **下游孤立 NB-ARC 段**：位点 B 只含 NB-ARC/linker/NA 类 motif，且落在某个「以 LRR
///    结尾」的位点 A 下游 `tail_gap` 之内——即串联结构中的尾随 NB-ARC（Prv 形态）。
///
/// 两类都只是把 B 的 motif 并入 A（`add_motifs` 会重新排序）并丢弃 B，不改变 A 的边界；
/// 基于 RefPlantNLR 400 条实测：位点 402 → 399，架构一致率 +3 例且零回退。
fn absorb_redundant_loci(
    nlrs: Vec<MotifList>,
    tail_gap: u64,
    def: &AnnotatorSignatureDefinition,
) -> Vec<MotifList> {
    let mut groups: HashMap<(String, Strand), Vec<usize>> = HashMap::new();
    let mut group_order: Vec<(String, Strand)> = Vec::new();
    for (i, list) in nlrs.iter().enumerate() {
        let key = (list.sequence_name().to_string(), list.strand());
        if !groups.contains_key(&key) {
            group_order.push(key.clone());
        }
        groups.entry(key).or_default().push(i);
    }

    let mut lists: Vec<Option<MotifList>> = nlrs.into_iter().map(Some).collect();
    for key in &group_order {
        let mut idx = groups[key].clone();
        idx.sort_by_key(|&i| lists[i].as_ref().map(|l| l.span().0).unwrap_or(0));
        let mut dropped = vec![false; idx.len()];

        // 规则 1：完全嵌套
        for a in 0..idx.len() {
            if dropped[a] {
                continue;
            }
            let (a_start, a_end) = lists[idx[a]].as_ref().unwrap().span();
            for b in (a + 1)..idx.len() {
                if dropped[b] {
                    continue;
                }
                let (b_start, b_end) = lists[idx[b]].as_ref().unwrap().span();
                if b_start >= a_end {
                    break; // 组内按起点排序，后续位点不可能再嵌套
                }
                if b_start >= a_start && b_end <= a_end {
                    let inner = lists[idx[b]].take().unwrap();
                    lists[idx[a]].as_mut().unwrap().add_motifs(&inner.motifs);
                    dropped[b] = true;
                }
            }
        }

        // 规则 2：下游孤立 NB-ARC 段（紧邻以 LRR 结尾的位点）
        for a in 0..idx.len() {
            if dropped[a] {
                continue;
            }
            if !lists[idx[a]]
                .as_ref()
                .unwrap()
                .domain_string(def)
                .ends_with("LRR")
            {
                continue;
            }
            let (a_start, a_end) = lists[idx[a]].as_ref().unwrap().span();
            let is_forward = lists[idx[a]].as_ref().unwrap().is_forward();
            let tail_candidates: Vec<(usize, u64)> = if is_forward {
                ((a + 1)..idx.len())
                    .filter(|&b| !dropped[b])
                    .filter_map(|b| {
                        let (b_start, _) = lists[idx[b]].as_ref().unwrap().span();
                        (b_start >= a_end).then_some((b, b_start - a_end))
                    })
                    .collect()
            } else {
                (0..a)
                    .rev()
                    .filter(|&b| !dropped[b])
                    .filter_map(|b| {
                        let (_, b_end) = lists[idx[b]].as_ref().unwrap().span();
                        (b_end <= a_start).then_some((b, a_start - b_end))
                    })
                    .collect()
            };
            for (b, tail_gap_bp) in tail_candidates {
                if dropped[b] {
                    continue;
                }
                let tail = lists[idx[b]].as_ref().unwrap();
                if tail_gap_bp > tail_gap {
                    break;
                }
                let only_nbarc = tail.motifs.iter().all(|m| {
                    matches!(
                        def.category(m.id),
                        DomainCategory::Nbarc | DomainCategory::Linker | DomainCategory::Na
                    )
                }) && tail.motifs.iter().any(|m| def.is_nbarc(m.id));
                if only_nbarc {
                    let seg = lists[idx[b]].take().unwrap();
                    lists[idx[a]].as_mut().unwrap().add_motifs(&seg.motifs);
                    dropped[b] = true;
                }
            }
        }
    }
    lists.into_iter().flatten().collect()
}

fn elongate(
    merged_seeds: Vec<MotifList>,
    motifs: &[Motif],
    distance: u64,
    def: &AnnotatorSignatureDefinition,
) -> Vec<MotifList> {
    // usedMotifs dedup key: `(id, dna_start)`.
    let mut used: HashSet<(u8, u64)> = HashSet::new();
    let mut pre_nlrs: Vec<MotifList> = Vec::new();

    for mut list in merged_seeds {
        // Entry check: skip if the seed contains an already-used motif.
        if list
            .motifs
            .iter()
            .any(|m| used.contains(&(m.id, m.dna_start)))
        {
            continue;
        }

        // Locate the index of the original firstMotif within the full list (equals = id + dnaStart).
        let first = list.first_motif().clone();
        let first_index = motifs
            .iter()
            .position(|m| m.id == first.id && m.dna_start == first.dna_start)
            .unwrap_or(0);

        // Elongate forward (index decreasing, firstMotif changes dynamically).
        let mut idx = first_index as isize - 1;
        while idx >= 0 {
            let motif = &motifs[idx as usize];
            let cur_first = list.first_motif().clone();
            if cur_first.dna_start.abs_diff(motif.dna_start) <= distance {
                if def.rank(motif.id) < def.rank(cur_first.id) {
                    let key = (motif.id, motif.dna_start);
                    if !used.contains(&key) && (!motif.has_stop() || def.is_nbarc(motif.id)) {
                        list.add_motif(motif.clone());
                    }
                }
            } else {
                break;
            }
            idx -= 1;
        }

        // Elongate backward (re-locate the dynamic lastMotif from first_index, index increasing).
        let mut idx = first_index;
        while idx < motifs.len() {
            let cur_last = list.last_motif().clone();
            if motifs[idx].id == cur_last.id && motifs[idx].dna_start == cur_last.dna_start {
                break;
            }
            idx += 1;
        }
        let mut idx = idx + 1;
        while idx < motifs.len() {
            let motif = &motifs[idx];
            let cur_last = list.last_motif().clone();
            if cur_last.dna_start.abs_diff(motif.dna_start) <= distance {
                let r1 = def.rank(motif.id);
                let r2 = def.rank(cur_last.id);
                if r1 > r2 || (r1 == r2 && def.is_lrr(motif.id)) {
                    let key = (motif.id, motif.dna_start);
                    if !used.contains(&key) && (!motif.has_stop() || def.is_nbarc(motif.id)) {
                        list.add_motif(motif.clone());
                    }
                }
            } else {
                break;
            }
            idx += 1;
        }

        // Mark as used.
        for m in &list.motifs {
            used.insert((m.id, m.dna_start));
        }
        pre_nlrs.push(list);
    }

    pre_nlrs
}
