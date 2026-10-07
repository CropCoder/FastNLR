# nlr-scan

> Sliding-window motif scanning: PWM scoring, p-value filtering and non-overlap arbitration.

`nlr-scan` answers one question per translated protein sequence: *which motif hits are
statistically significant, and which of the overlapping candidates wins?* It slides every
motif's position-weight matrix across the sequence, converts raw scores into p-values through
the CDF table, and arbitrates between overlapping hits so that the returned list contains a
non-overlapping set of the most significant motifs.

The hot loop is designed around an integer pre-filter and 8-lane SIMD, so the great majority of
windows never reach a floating-point lookup.

## Position in the workspace

| | |
|---|---|
| Layer | 3 |
| Depends on | `nlr-core`, `nlr-config`, `nlr-seq`, `wide` |
| Used by | `nlr-cli` |
| Input | a translated protein sequence (`&str`) |
| Output | `MotifList` with protein-side hits |

## Algorithm

1. **Iterate motifs in profile order.** Motifs are processed in the order they first appear in
   `mot.txt`, which keeps results independent of hash-map iteration order.
2. **Score each window.** The PWM contributes integer scores, summed over the motif width.
3. **Pre-filter cheaply.** `MotifDefinition::score_thresholds` yields an integer threshold per
   motif; windows below it cannot reach the preliminary p-value and are discarded without a
   floating-point operation.
4. **Convert to a p-value.** Surviving windows are looked up in the right-tailed CDF.
5. **Arbitrate overlaps.** A candidate is compared against already-accepted hits whose spans
   intersect it, using the product of their p-values; the weaker candidate is displaced and the
   stronger one is kept.
6. **Apply the acceptance threshold.** Only hits with `p < accept` are materialised as `Motif`
   records.

Steps 2–3 run eight windows at a time via `wide::i32x8`, with a scalar tail loop for the
remaining windows.

### Thresholds

| Constant | Default | Role |
|---|---|---|
| `THRESH_PRELIMINARY` | `1e-4` | integer pre-filter cut-off — cheap rejection of clearly insignificant windows |
| `THRESH_ACCEPT` | `1e-5` | final acceptance cut-off for a motif hit |

Both are overridable through `MotifParser::with_thresholds`, which `fastnlr` exposes as
`--motif-prelim-p` and `--motif-accept-p`. Raising the acceptance threshold (for example to
`1e-3`) increases recall for highly divergent helper NLRs at the cost of more weak hits.

## Public API

| Item | Description |
|---|---|
| `MotifParser::new(definition)` | parser with the default thresholds |
| `MotifParser::with_thresholds(definition, prelim, accept)` | parser with custom thresholds |
| `MotifParser::find_motifs(protein_seq_id, protein)` | scans one protein sequence and returns its accepted hits |
| `has_signature(list, def)` | convenience wrapper over `AnnotatorSignatureDefinition::has_signature` |

`find_motifs` returns an empty `MotifList` when the sequence is shorter than the longest motif,
and materialises each accepted hit with:

- `position` — 1-based offset of the hit in the protein sequence;
- `protein_sequence` — the matched substring, i.e. the motif width;
- `pvalue` — the arbitrated (non-overlap) p-value;
- `protein_sequence_id` — propagated from the caller, later parsed by `nlr-cli` to recover the
  contig, fragment offset and reading frame.

## Example

```rust
use nlr_config::MotifDefinition;
use nlr_scan::MotifParser;

let def = MotifDefinition::load_from_str(mot_txt, store_txt).unwrap();
let parser = MotifParser::new(def);

let hits = parser.find_motifs("chr1_0_frame+0", &protein);
for m in &hits.motifs {
    println!("motif_{} at {} p={:e}", m.id, m.position, m.pvalue);
}
```

## Guarantees

- Output order is deterministic: motifs are scanned in profile order and the arbitration is
  reproducible for identical input.
- Overlapping hits are resolved, never returned together, so downstream assembly never has to
  deal with duplicate coverage of the same residues.
- A window is scored with the full motif width; partial windows at the sequence end are never
  scored.

## Testing

```bash
cargo test -p nlr-scan
```

`tests/scan.rs` covers score accumulation, threshold behaviour, overlap arbitration and the
SIMD/scalar boundary.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
