# nlr-report

> Run statistics aggregation for FastNLR.

`nlr-report` condenses a finished run into a small, deterministic statistics model and renders
it as TSV. It answers the questions users ask first — how many motifs and loci were found,
how many are complete, how the counts distribute over contigs, and which motifs dominate.

The crate has no I/O of its own: it aggregates data structures and writes to any
`std::io::Write`.

## Position in the workspace

| | |
|---|---|
| Layer | 4 |
| Depends on | `nlr-core` |
| Used by | `nlr-cli` (`.stats.tsv`), `nlr-plot` (chart input) |
| Input | motif hits grouped by contig, assembled loci, signature definition |
| Output | `RunStats` and a TSV rendering of it |

## Public API

| Item | Description |
|---|---|
| `collect(motifs_by_seq, nlrs, def) -> RunStats` | aggregates global, per-contig, per-motif and per-type counts |
| `write_tsv(w, stats) -> io::Result<()>` | renders the statistics as a sectioned TSV report |

### `RunStats`

| Field | Meaning |
|---|---|
| `seq_count` | number of contigs that produced at least one motif |
| `motif_total` | total accepted motif hits |
| `nlr_total` | total assembled loci |
| `nlr_complete` | loci satisfying `is_complete_nlr` (P-loop plus at least one LRR) |
| `per_chromosome` | `contig -> (motifs, loci, complete loci)` |
| `motif_counts` | `motif id -> hit count` |
| `nlr_type_counts` | `domain string -> locus count`, e.g. `CC-NBARC-LRR` |

All maps are `BTreeMap`, so iteration order — and therefore the output order — is deterministic
regardless of the order in which contigs were processed in parallel.

## TSV layout

```
#section	key	value
global	sequence_count	12
global	motif_total	2403
global	nlr_total	57
global	nlr_complete	41
chromosome	chr1	812	21	17
motif	motif_1	61
```

| Section | Columns |
|---|---|
| `global` | `key`, `value` |
| `chromosome` | contig name, motif count, locus count, complete locus count (extends the 3-column header to 5 fields) |
| `motif` | `motif_<id>`, hit count |

## Notes

- The header declares three columns while `chromosome` rows carry five; consumers should key off
  the section name rather than the column count.
- `nlr_type_counts` is populated for `nlr-plot` and is not emitted by `write_tsv` in this
  version.

## Example

```rust
use nlr_report::{collect, write_tsv};

let stats = collect(&result.motifs_by_seq, &result.nlrs, &result.def);
println!("{} motifs, {} loci ({} complete)", stats.motif_total, stats.nlr_total, stats.nlr_complete);

let mut buf = Vec::new();
write_tsv(&mut buf, &stats).unwrap();
```

## Testing

```bash
cargo test -p nlr-report
```

The crate currently has no dedicated test suite; it is exercised indirectly through the
`nlr-cli` pipeline tests.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
