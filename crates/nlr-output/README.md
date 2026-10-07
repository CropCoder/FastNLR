# nlr-output

> Multi-format writers for FastNLR results.

`nlr-output` renders assembled loci into the formats downstream users actually consume: a plain
text report, GFF, BED12, motif-level BED, an NB-ARC multiple alignment, extracted locus
sequences and a TSV motif exchange format. Every writer is generic over `std::io::Write`, so the
same code produces files, in-memory buffers or test fixtures.

The crate's central concern is **coordinate consistency**: the intervals written to GFF, BED and
the extracted locus FASTA must describe exactly the same bases.

## Position in the workspace

| | |
|---|---|
| Layer | 4 |
| Depends on | `nlr-core` |
| Used by | `nlr-cli` |
| I/O | writes to any `Write`; the locus extractor receives contig sequences from the caller |

## Writers

| Function | Output | Notes |
|---|---|---|
| `write_report_txt` | `.nlr.txt` | headerless, 7 columns: `seqname`, `name`, domain class, `start`, `end`, strand, motif list |
| `write_nlr_gff` | `.nlr.gff` | 5 comment lines, feature `NBSLRR`, source `FastNLR`, `start + 1` (GFF is 1-based), `name=` and `nlrClass=` attributes; optional `complete_only` filter |
| `write_nlr_bed` | `.nlr.bed` | BED12 with `itemRgb`; one block per motif; blocks are emitted in reverse order on the reverse strand |
| `write_motif_bed` | `.motifs.bed` | 9 columns per motif, p-value in the score column, per-motif colour; optional `STOP` annotation rows |
| `write_nbarc_alignment_fasta` | `.nbarc.fasta` | P-loop-anchored alignment, gap-padded per missing consensus block, `*`/`X` replaced by `_`; optional CED-4 reference row |
| `write_nlr_loci` | locus FASTA | extraction from a single contig (legacy entry point) |
| `write_nlr_loci_all` | `.loci.fasta` | extraction across all contigs, ±`--flank` bases, reverse-complemented on the reverse strand, wrapped at 100 bp |
| `export_motifs` / `import_motifs` | TSV | 11-field round trip used by the checkpoint mechanism |

### BED colour convention

| Colour | Meaning |
|---|---|
| `0,255,0` green | complete NLR (P-loop and at least one LRR) |
| `255,128,0` orange | partial NLR |
| `255,0,0` red | the locus contains a stop codon (takes priority over completeness) |

### Locus FASTA header

```
>{locus_name} {contig} {start}-{end} strand:{+|-} {motif_list}
```

The coordinates in the header are the same clamped interval used to slice the sequence, so the
extracted bases can be reproduced directly from the GFF/BED interval.

## Coordinate guarantees

1. The alignment FASTA locates the P-loop with the corrected `while(!is_ploop)` rule.
2. The alignment FASTA replaces stop codons and unknown amino acids with `_`.
3. Locus extraction applies a unified clamp at the end of every contig (`saturating_sub` on the
   left, `min(contig_length)` on the right) instead of assuming enough flanking sequence.
4. Locus extraction iterates **all** contigs, not only the first.
5. Extracted sequence coordinates always match the reported GFF/BED coordinates.
6. Small p-values are written in scientific notation rather than being flushed to zero.

## Example

```rust
use nlr_output::{write_nlr_gff, write_nlr_loci_all};

let mut gff = Vec::new();
write_nlr_gff(&mut gff, &loci, &def, "2026-10-07 12:00:00", false).unwrap();

let contigs: Vec<(&str, &str)> = vec![("chr1", chr1_sequence)];
let mut fasta = Vec::new();
write_nlr_loci_all(&mut fasta, &loci, &contigs, 2_000).unwrap();
```

## Notes

- Writers never open files; path handling, output naming and directory creation belong to
  `nlr-cli`.
- `write_nlr_loci` (single contig) is retained for API compatibility; `fastnlr` uses
  `write_nlr_loci_all`.
- `export_motifs` / `import_motifs` are the serialisation format behind
  `--checkpoint`, and are the only place where a locus-free motif list round-trips through text.

## Testing

```bash
cargo test -p nlr-output
```

`tests/output.rs` covers every writer, including header layout, block ordering on the reverse
strand, colour selection and the clamp behaviour at contig ends.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
