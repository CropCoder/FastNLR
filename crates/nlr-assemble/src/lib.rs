//! nlr-assemble — three-step NLR locus assembly: findSeeds -> mergeSeeds -> elongate.
//!
//! Faithful reimplementation of Java `NLR_Annotator.findNLRs` (activated version only):
//! - findSeeds: build seed combinations from adjacent motifs within
//!   `distanceWithinMotifCombination`;
//! - mergeSeeds: same strand + 4 dna_start combination differences within threshold + rank
//!   monotonic;
//! - elongate: extend at both ends by rank rules, with usedMotifs dedup.

use std::collections::HashSet;
use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;
use nlr_core::signature_def::AnnotatorSignatureDefinition;

/// Three-step assembly parameters.
#[derive(Debug, Clone, Copy)]
pub struct AssembleParams {
    /// Maximum distance between adjacent motifs within a seed (default 500).
    pub distance_within_motif_combination: u64,
    /// Elongation distance (default 2500).
    pub distance_for_elongating: u64,
    /// Seed merge distance (default 10000).
    pub distance_between_motif_combinations: u64,
}

impl Default for AssembleParams {
    fn default() -> Self {
        AssembleParams {
            distance_within_motif_combination: 500,
            distance_for_elongating: 2500,
            distance_between_motif_combinations: 10000,
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

    let seeds = find_seeds(&motifs, params.distance_within_motif_combination, def);
    let merged = merge_seeds(seeds, params.distance_between_motif_combinations, def);
    let mut pre_nlrs = elongate(merged, &motifs, params.distance_for_elongating, def);

    // Java sorts `preNlrs` before naming (forward strand ascending / reverse strand
    // descending), which determines the nlrN naming order.
    pre_nlrs.sort();

    let mut nlrs: Vec<MotifList> = Vec::new();
    let mut suffix = 1usize;
    for mut list in pre_nlrs {
        list.sort();
        list.name = format!("{}_nlr{}", seq_id, suffix);
        suffix += 1;
        nlrs.push(list);
    }
    nlrs
}

/// (1) findSeeds (Java `findNLRs_substep1_findSeeds`).
fn find_seeds(
    motifs: &[Motif],
    distance: u64,
    def: &AnnotatorSignatureDefinition,
) -> Vec<MotifList> {
    let mut seeds = Vec::new();
    let mut count = 0usize;

    let mut i = 0;
    while i < motifs.len() {
        let mut potential = vec![motifs[i].clone()];
        let mut s: Vec<u8> = vec![motifs[i].id];

        let mut j = i + 1;
        let mut cur = &motifs[i];
        while j < motifs.len() {
            let next = &motifs[j];
            if cur.dna_start.abs_diff(next.dna_start) <= distance {
                s.push(next.id);
                potential.push(next.clone());
                cur = next;
                if def.is_seed(&s) {
                    count += 1;
                    let name = format!("{}_nlr{}", motifs[i].dna_sequence_id.as_deref().unwrap_or(""), count);
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

/// (2) mergeSeeds (Java `findNLRs_subsep2a_mergeSeeds`).
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

/// (3) elongate (Java `findNLRs_substep3_elongate`, activated version only).
///
/// Key replication point: inside the loop Java uses `motifList.firstMotif()` /
/// `motifList.lastMotif()` as a **dynamic** baseline (addMotif re-sorts, so the baseline
/// shifts as elongation proceeds). The first/last from the initial clone must not be pinned.
fn elongate(
    merged_seeds: Vec<MotifList>,
    motifs: &[Motif],
    distance: u64,
    def: &AnnotatorSignatureDefinition,
) -> Vec<MotifList> {
    // usedMotifs dedup key: `(id, dna_start)`, equivalent to Java `Motif.equals`.
    let mut used: HashSet<(u8, u64)> = HashSet::new();
    let mut pre_nlrs: Vec<MotifList> = Vec::new();

    for mut list in merged_seeds {
        // Entry check: skip if the seed contains an already-used motif.
        if list.motifs.iter().any(|m| used.contains(&(m.id, m.dna_start))) {
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
