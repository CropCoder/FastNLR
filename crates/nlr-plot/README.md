# nlr-plot

> PNG statistics charts for FastNLR runs.

`nlr-plot` renders the run statistics produced by `nlr-report` into three bitmaps:
motif hit counts, per-contig locus counts and locus type counts. It uses `plotters` with the
`ab_glyph` font backend, which has no system font dependency, so the crate compiles statically
and cross-compiles cleanly.

## Position in the workspace

| | |
|---|---|
| Layer | 5 |
| Depends on | `nlr-report`, `plotters` |
| Used by | `nlr-cli` (`.plots/` directory) |
| Output | PNG files, 1200 × 800 px |

## Charts

| Function | File written by `fastnlr` | Contents |
|---|---|---|
| `plot_motif_counts` | `01-motif-counts.png` | bar chart of hits per motif id; the x-axis range is derived from the data, so the 20-motif built-in profile and the 28-motif extended library both render correctly |
| `plot_chromosome_nlrs` | `02-chromosome-nlrs.png` | bar chart of loci per contig, labelled with contig names |
| `plot_nlr_type_counts` | `03-nlr-types.png` | bar chart of locus types (for example `CC-NBARC-LRR`), sorted by count descending |

All three functions have the same signature shape:

```rust
pub fn plot_...(path: &Path, stats: &RunStats) -> Result<(), Box<dyn std::error::Error>>
```

and take their data from `nlr_report::RunStats`.

## Font handling

The `ab_glyph` backend used by `plotters` ships **no** fonts, and rendering text without one
fails with `FontUnavailable` — the historical cause of blank white PNGs. This crate therefore
embeds `assets/LiberationSans-Regular.ttf` (SIL Open Font License 1.1, compatible with this
project's GPL-3.0 licensing) and registers it once, guarded by a `std::sync::Once`, as the
`sans-serif` family.

No system font configuration is read at build time or at run time.

## Rendering notes

- The x-axis of the motif and contig charts uses a segmented coordinate so that bars are centred
  on their category rather than on a numeric position.
- Contig labels are formatted through an axis formatter that handles `Exact`, `CenterOf` and
  `Last` segment values; the number of printed labels is capped to avoid overlap on assemblies
  with many contigs.
- Charts are written as PNG. The bitmap backend is intentional: it keeps the dependency tree
  free of FreeType/fontconfig, which matters for static and cross-compiled builds.

## Example

```rust
use nlr_plot::{plot_motif_counts, plot_nlr_type_counts};

plot_motif_counts(Path::new("01-motif-counts.png"), &stats).unwrap();
plot_nlr_type_counts(Path::new("03-nlr-types.png"), &stats).unwrap();
```

## Testing

```bash
cargo test -p nlr-plot
```

The crate currently has no dedicated test suite. Rendering was validated manually by asserting
that the produced PNGs are non-blank and mutually different; adding an automated check is
tracked as future work.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

The bundled Liberation Sans font is licensed under the SIL Open Font License 1.1.

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
