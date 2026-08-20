//! nlr-scan — scanning engine: sliding-window scoring + non-overlap arbitration + signature pre-filtering.
//!
//! Faithful reimplementation of Java `MotifParser.findMotifs` and `get_motifs_stats`:
//! - For each protein sequence, scans motifs in order of first appearance;
//! - Scoring: PWM integer accumulation + `(int)` truncation + CDF table lookup;
//! - Preliminary threshold `p < 1e-4`, final acceptance `p < 1e-5`;
//! - Non-overlap arbitration: comparison via product of p-values.
//!
//! Performance optimizations:
//! - **Score threshold table**: precompute `kmer_score < T[id] ⟺ p < 1e-4`; the vast majority of
//!   windows are filtered by pure integer comparison, and only windows passing the preliminary
//!   filter query the CDF for an exact p-value;
//! - **SIMD cross-window scoring**: use `wide` to accumulate integer scores for 8 consecutive
//!   windows in parallel.

use nlr_config::MotifDefinition;
use nlr_core::motif::Motif;
use nlr_core::motif_list::MotifList;

/// Preliminary threshold (Java `thresh = 0.0001D`).
pub const THRESH_PRELIMINARY: f64 = 1e-4;
/// Final acceptance threshold (Java `1E-5`).
pub const THRESH_ACCEPT: f64 = 1e-5;

/// Number of SIMD parallel windows (wide i32x8).
const SIMD_LANES: usize = 8;

/// Scanning engine.
pub struct MotifParser {
    definition: MotifDefinition,
    /// Preliminary integer threshold: when `score < prelim_thresholds[id]`, `p < 1e-4`.
    prelim_thresholds: Vec<i32>,
    /// Final integer threshold: when `score < accept_thresholds[id]`, `p < 1e-5`.
    /// (Currently find_motifs judges final acceptance via exact pvalue; this table is
    /// retained for a future fast path.)
    #[allow(dead_code)]
    accept_thresholds: Vec<i32>,
}

impl MotifParser {
    pub fn new(definition: MotifDefinition) -> Self {
        let prelim_thresholds = definition.score_thresholds(THRESH_PRELIMINARY);
        let accept_thresholds = definition.score_thresholds(THRESH_ACCEPT);
        MotifParser {
            definition,
            prelim_thresholds,
            accept_thresholds,
        }
    }

    /// Scan all motifs against a protein sequence (equivalent to `findMotifs`).
    ///
    /// Returns the MotifList for the sequence; returns empty when the sequence is shorter
    /// than the maximum motif length.
    pub fn find_motifs(&self, protein_seq_id: &str, protein: &str) -> MotifList {
        let mut motif_list = MotifList::new(protein_seq_id.to_string(), Vec::new());

        let seq_len = protein.len();
        if seq_len < self.definition.max_length() as usize {
            return motif_list;
        }

        let (hits, pvalues) = self.get_motifs_stats(protein, seq_len);

        let bytes = protein.as_bytes();
        for (q, &hit) in hits.iter().enumerate() {
            if hit == 0 {
                continue;
            }
            let motif_id = hit; // hit directly stores the motif id (1..=20)
            let p = pvalues[q];
            if p < THRESH_ACCEPT {
                let width = self.definition.length(motif_id) as usize;
                let seq = std::str::from_utf8(&bytes[q..q + width])
                    .unwrap_or("")
                    .to_string();
                let m = Motif::new_protein(
                    motif_id,
                    protein_seq_id.to_string(),
                    q as u64 + 1, // 1-based
                    seq,
                    p,
                );
                motif_list.add_motif(m);
            }
        }

        motif_list
    }

    /// Non-overlap arbitration (equivalent to `get_motifs_stats`); returns `(hits, pvalues)`.
    ///
    /// - `hits[q]` = the motif id hit at this position (0 = none);
    /// - `pvalues[q]` = the non-overlap p-value at this position (initial 1.0).
    fn get_motifs_stats(&self, protein: &str, seq_len: usize) -> (Vec<u8>, Vec<f64>) {
        let mut hits = vec![0u8; seq_len];
        let mut pvalues = vec![1.0f64; seq_len];
        let mut maxws: usize = 0;
        let bytes = protein.as_bytes();

        for &motif_id in self.definition.motif_names() {
            let ws = self.definition.length(motif_id) as usize;
            maxws = maxws.max(ws);

            if seq_len < ws {
                continue;
            }

            let prelim_t = self.prelim_thresholds[motif_id as usize];
            let j_end = seq_len - ws + 1;

            // SIMD batch scoring: compute scores for 8 windows at a time, pre-filter by pure
            // integer threshold. CDF right tail is decreasing: score >= prelim_t is significant
            // (p < 1e-4).
            let mut j = 0;
            while j + SIMD_LANES <= j_end {
                let scores = self.score_batch_simd(bytes, motif_id, ws, j);
                for lane in 0..SIMD_LANES {
                    let kmer_score = scores[lane];
                    if kmer_score < prelim_t {
                        // p >= 1e-4, not significant, skip.
                        continue;
                    }
                    let pos = j + lane;
                    let pvalue = self.definition.cdf(motif_id, kmer_score);
                    self.arbitrate(motif_id, pos, ws, pvalue, &mut hits, &mut pvalues, maxws, seq_len);
                }
                j += SIMD_LANES;
            }
            // Scalar handling for the remaining windows.
            while j < j_end {
                let kmer_score = self.score_at(bytes, motif_id, ws, j);
                if kmer_score >= prelim_t {
                    let pvalue = self.definition.cdf(motif_id, kmer_score);
                    self.arbitrate(motif_id, j, ws, pvalue, &mut hits, &mut pvalues, maxws, seq_len);
                }
                j += 1;
            }
        }

        (hits, pvalues)
    }

    /// Non-overlap arbitration (single window, extracted for reuse).
    #[inline]
    fn arbitrate(
        &self,
        motif_id: u8,
        j: usize,
        ws: usize,
        pvalue: f64,
        hits: &mut [u8],
        pvalues: &mut [f64],
        maxws: usize,
        seq_len: usize,
    ) {
        let mut ok_to_mark = true;
        let mut prod = 1.0f64;

        let first = j.saturating_sub(maxws - 1);
        let last = (j + ws).min(seq_len);

        // Look to the right.
        for k in j..last {
            if !ok_to_mark {
                break;
            }
            if hits[k] != 0 {
                prod *= pvalues[k];
                if pvalue >= prod {
                    ok_to_mark = false;
                }
            }
        }
        // Look to the left.
        for k in first..j {
            if !ok_to_mark {
                break;
            }
            let h = hits[k];
            if h != 0 && self.definition.length(h) as usize > j - k {
                prod *= pvalues[k];
                if pvalue >= prod {
                    ok_to_mark = false;
                }
            }
        }

        if ok_to_mark {
            hits[j] = motif_id;
            pvalues[j] = pvalue;
            for k in (j + 1)..last {
                hits[k] = 0;
                pvalues[k] = 1.0;
            }
            for k in first..j {
                let h = hits[k];
                if h != 0 && self.definition.length(h) as usize > j - k {
                    hits[k] = 0;
                    pvalues[k] = 1.0;
                }
            }
        }
    }

    /// Single-window scoring (scalar).
    #[inline]
    fn score_at(&self, bytes: &[u8], motif_id: u8, ws: usize, pos: usize) -> i32 {
        let mut score: i32 = 0;
        for i in 0..ws {
            score += self.definition.score(motif_id, i, bytes[pos + i]);
        }
        score
    }

    /// SIMD batch scoring: accumulate integer scores in parallel for 8 consecutive windows
    /// `[pos, pos+8)`.
    #[inline]
    fn score_batch_simd(&self, bytes: &[u8], motif_id: u8, ws: usize, pos: usize) -> [i32; SIMD_LANES] {
        use wide::i32x8;
        let mut acc = i32x8::splat(0);
        for i in 0..ws {
            let lane: [i32; SIMD_LANES] = std::array::from_fn(|l| {
                self.definition.score(motif_id, i, bytes[pos + l + i]) as i32
            });
            acc += i32x8::from(lane);
        }
        acc.to_array()
    }

    /// Single-window scoring + CDF (retained for paths that need an exact p-value).
    #[inline]
    #[allow(dead_code)]
    fn pval_at(&self, bytes: &[u8], motif_id: u8, ws: usize, pos: usize) -> f64 {
        let score = self.score_at(bytes, motif_id, ws, pos);
        self.definition.cdf(motif_id, score)
    }
}

/// Whether the NLR signature is present (equivalent to `MotifList.has_NLR_signature`).
pub fn has_signature(
    list: &MotifList,
    def: &nlr_core::signature_def::AnnotatorSignatureDefinition,
) -> bool {
    let ids: Vec<u8> = list.motifs.iter().map(|m| m.id).collect();
    def.has_signature(&ids)
}
