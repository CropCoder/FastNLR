# fastnlr (crate `nlr-cli`)

> Command-line entry point and pipeline orchestration for FastNLR.

This crate owns everything that is not an algorithm: argument parsing, input validation, the
parallel scan pipeline, checkpointing, graceful interruption, run summaries and output file
management. It links the seven library crates together and produces the `fastnlr` binary.

The package is published as `fastnlr`; the library target keeps the internal name `nlr_cli` so
that in-tree imports and integration tests stay stable.

## Position in the workspace

| | |
|---|---|
| Layer | 6 (top) |
| Depends on | `nlr-core`, `nlr-config`, `nlr-seq`, `nlr-scan`, `nlr-assemble`, `nlr-output`, `nlr-report`, `nlr-plot`, `clap`, `rayon`, `tracing`, `indicatif`, `ctrlc`, `tempfile` |
| Targets | binary `fastnlr`, library `nlr_cli` |
| Embeds | `data/mot.txt` and `data/store.txt` via `include_str!` |

## Usage

```bash
fastnlr -i genome.fasta -o results
fastnlr -i genome.fasta -o results -t 8
fastnlr -i genome.fasta -o results --checkpoint ckpt/ --flexible-seed 3
```

`-i` and `-o` are required; the built-in motif profiles are used unless `-x` / `-y` are given.

### Flag groups

| Group | Flags |
|---|---|
| Input | `-i <INPUT>` |
| Motif config | `-x <mot.txt>`, `-y <store.txt>` |
| Output | `-o <DIR>`, `--flank <bp>` (default 2000) |
| Performance | `-t <threads>`, `-n <fragments per chunk>` (default 1000) |
| Run control | `--checkpoint <DIR>`, `--tmpdir <DIR>` |
| Observability | `--progress <auto\|bar\|simple\|off>`, `--log-level <level>` |
| Recall | `--motif-accept-p`, `--motif-prelim-p`, `--relaxed-seed <n>`, `--flexible-seed <n>`, `--require-ploop` |
| Motif library metadata | `--motif-category <ID=CAT>`, `--seed-combination <IDS>`, `--signature <IDS>` |

## Pipeline

`run` / `run_with_progress` execute the following stages:

1. **Load the profile and build the rule definition.** User paths take precedence over the
   embedded defaults; CLI-declared categories, seeds, signatures and the flexible-seeding
   configuration are injected into `AnnotatorSignatureDefinition`.
2. **Checkpoint short-circuit.** If `--checkpoint` names a directory containing `motifs.tsv`,
   the scan is skipped and only assembly runs.
3. **Stream and scan.** `SequenceChopper` yields fragments in fixed-size chunks; each chunk is
   scanned in parallel by a dedicated rayon pool with one task per fragment. The chunk is
   dropped after merging, so peak memory during this stage is bounded by the chunk size rather
   than by the genome.
4. **Normalise and pre-filter.** Per contig, hits are sorted and deduplicated, then filtered by
   local motif clusters grouped by `(contig, strand, frame)`. Cluster boundaries come from motif
   coordinates, not from fragment boundaries, which makes results independent of the chopping
   phase.
5. **Assemble.** Contigs are assembled in parallel and re-ordered by sorted contig id so that
   output order is deterministic.
6. **Write the checkpoint** (if requested) and return a `RunResult`.

`SIGINT` sets an atomic flag; the scan stops after the current chunk and the completed results
are still written out.

## Public API

| Item | Description |
|---|---|
| `RunConfig` | full run configuration, see below |
| `RunResult` | `nlrs`, `motifs_by_seq`, `def` |
| `run(&RunConfig)` / `run_with_progress(&RunConfig, Option<&ProgressBar>)` | pipeline entry points; the second reports progress |
| `all_motifs(&RunResult)` | flattens motif hits across contigs |
| `save_checkpoint(dir, &RunResult)` / `load_checkpoint(dir, def)` | TSV-based checkpoint round trip |
| `EMBEDDED_MOT` / `EMBEDDED_STORE` | the built-in profiles |

### `RunConfig`

| Field | Default | Source flag |
|---|---|---|
| `input_fasta` | — | `-i` |
| `mot_file` / `store_file` | `None` (embedded) | `-x` / `-y` |
| `fragment_length` | 20 000 | internal |
| `overlap` | 2 000 | internal |
| `threads` | available parallelism | `-t` |
| `seqs_per_thread` | 1 000 | `-n` |
| `checkpoint_dir` | `None` | `--checkpoint` |
| `assemble` | `AssembleParams::default()` | internal |
| `motif_accept_p` | `1e-5` | `--motif-accept-p` |
| `motif_prelim_p` | `1e-4` | `--motif-prelim-p` |
| `motif_categories` | empty | `--motif-category` |
| `extra_seeds` / `extra_signatures` | empty | `--seed-combination` / `--signature` |
| `flexible_seed` | `None` | `--flexible-seed` / `--require-ploop` |

## Outputs

For an input `genome.fasta`, the prefix `genome` is derived by stripping `.fa/.fasta/.fna/.fas`
and `.gz`, and the following files are written into `-o`:

| File | Writer |
|---|---|
| `genome.nlr.txt` | `nlr-output` |
| `genome.nlr.gff` | `nlr-output` |
| `genome.nlr.bed` | `nlr-output` |
| `genome.motifs.bed` | `nlr-output` |
| `genome.nbarc.fasta` | `nlr-output` |
| `genome.loci.fasta` | `nlr-output` |
| `genome.stats.tsv` | `nlr-report` |
| `genome.summary.txt` | this crate |
| `genome.plots/{01-motif-counts,02-chromosome-nlrs,03-nlr-types}.png` | `nlr-plot` |

All outputs are written by default; there is no opt-in flag per format.

## Notes

- `load_checkpoint` rebuilds the pipeline from the motif TSV and assembles with
  `AssembleParams::default()`, so assembly parameters set through `RunConfig` are not replayed
  on a checkpoint-only run.
- Locus extraction reads the genome through `nlr_seq::fasta::read_all`, so peak memory for a
  full run scales with the genome size during the output stage.
- The GFF date header uses UTC.

## Testing

```bash
cargo test -p fastnlr
```

`tests/run.rs` drives the library pipeline from synthetic DNA and asserts that motifs are
recovered; `tests/e2e.rs` covers the scan-to-assembly path on back-translated constructs. The
argument parser and the output-file set are not covered by automated tests yet.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
