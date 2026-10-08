# nlr-config

> Parser for the FastNLR motif profiles: `mot.txt` (PWM) and `store.txt` (CDF).

`nlr-config` turns the two plain-text motif profile files into a single in-memory
`MotifDefinition` that the scanning engine can query in constant time: an integer
position-weight matrix, a cumulative-distribution table for p-values, per-motif lengths and
precomputed score thresholds.

The crate is deliberately size-agnostic: the number of motifs is derived from the largest id
found in the files, so the built-in 20-motif profile and the 28-motif extended library
(`data/motif_library_rnl_tir/`) are handled by the same code path.

## Position in the workspace

| | |
|---|---|
| Layer | 2 |
| Depends on | `nlr-core`, `memmap2` |
| Used by | `nlr-scan`, `nlr-cli` |
| I/O | reads files (read-only, memory-mapped) or parses in-memory strings |

## Input formats

### `mot.txt` — position weight matrix

```
motif_1@0@G 12
motif_1@1@K 8
motif_1@2@T 5
```

One entry per line: `motif@position@amino_acid score`, where

- `motif` is `motif_<id>`, id starting at 1;
- `position` is **0-based**; motif length is `max(position) + 1`;
- `amino_acid` is a single upper-case letter;
- `score` is an integer.

The matrix is stored flattened as `[aa_index * width + position]` with `aa_index = aa - b'A'`,
so the per-entry stride is the motif width.

### `store.txt` — cumulative distribution function

```
motif_1@0 1.0
motif_1@30 1e-7
```

One entry per line: `motif@score pvalue`. The table is a **right-tail** CDF:
`p = P(random score >= score)`, monotonically decreasing in score.

Two under-specified cases are handled conservatively, so that an incomplete custom table can
never fabricate significant hits:

- scores **beyond the table's upper bound** return the table's last (smallest) p-value, which is
  a valid upper bound for the true p-value;
- score slots **missing from the table** default to `p = 1.0` (not significant) instead of `0.0`.

In addition, `load` / `load_from_str` reject a profile pair whose files disagree on the motif
ids, rather than panicking later during scoring.

## Public API

| Item | Description |
|---|---|
| `MotifDefinition::load(pwm_file, cdf_file)` | memory-maps and parses both files |
| `MotifDefinition::load_from_str(mot, store)` | parses in-memory text; this is how `nlr-cli` embeds the default profiles with `include_str!` |
| `load_default(dir)` | convenience wrapper for `<dir>/mot.txt` + `<dir>/store.txt` |
| `max_motif_id()` | highest id present in `mot.txt`; internal arrays are sized `max_motif_id + 1` |
| `motif_names()` | ids in order of first appearance — this is the scan iteration order |
| `length(id)` / `max_length()` | motif width in amino acids |
| `score(id, position, aa)` | PWM lookup; returns `0` for characters below `A` or outside `A..Z` |
| `cdf(id, score)` | p-value lookup; out-of-range scores return the last tabulated value |
| `score_thresholds(thresh)` | per-motif integer score `T` such that `score >= T` implies `p < thresh`; `i32::MAX` means no score can reach the threshold |
| `parse_motif_line(line)` | parses one exported TSV motif row into `(id, protein id, position, pvalue)` |

`MotifId` is re-exported as a `u8` alias.

## Performance notes

`score_thresholds` exists so that the scanning hot loop can reject the overwhelming majority of
windows with a single integer comparison, before any floating-point CDF lookup. Thresholds are
computed once per run, at parser construction.

File-backed loading uses `memmap2`, so a large profile is never copied into the heap; it is
parsed directly out of the mapped region.

## Example

```rust
use nlr_config::MotifDefinition;

let mot = "motif_1@0@G 12\nmotif_1@1@K 8\n";
let store = "motif_1@0 1.0\nmotif_1@20 1e-9\n";

let def = MotifDefinition::load_from_str(mot, store).unwrap();
assert_eq!(def.length(1), 2);
assert_eq!(def.score(1, 0, b'G'), 12);
assert!(def.cdf(1, 20) < 1e-5);
```

## Guarantees

- Parsing never panics on malformed input: unknown motif names, short lines and unparsable
  numeric fields are skipped.
- Ids may be sparse; absent ids become zero-length motifs and are effectively inert.
- The CDF is right-tailed and decreasing — callers must not treat it as a left-tail probability.
- `mot.txt` and `store.txt` must cover the same motif ids; a motif declared without a CDF table
  is reported as an `InvalidData` error.

## Testing

```bash
cargo test -p nlr-config
```

`tests/config.rs` covers multi-motif profiles, ids above 20, sparse ids and threshold
computation.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
