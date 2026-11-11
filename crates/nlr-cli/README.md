# fastnlr (crate `nlr-cli`)

> Command-line entry point and pipeline orchestration for FastNLR.

This crate owns everything that is not an algorithm: argument parsing, input validation, the
parallel scan pipeline, checkpointing, graceful interruption, run summaries and output file
management. It links the eight library crates together and produces the `fastnlr` binary.

The package is published as `fastnlr`; the library target keeps the internal name `nlr_cli` so
that in-tree imports and integration tests stay stable.

## Position in the workspace

| | |
|---|---|
| Layer | 6 (top) |
| Depends on | `nlr-core`, `nlr-config`, `nlr-seq`, `nlr-scan`, `nlr-assemble`, `nlr-domain`, `nlr-output`, `nlr-report`, `clap`, `rayon`, `tracing`, `indicatif`, `ctrlc`, `tempfile` |
| Targets | binary `fastnlr`, library `nlr_cli` |
| Embeds | `data/mot.txt` and `data/store.txt` via `include_str!` |

## Usage

```bash
fastnlr -i genome.fasta -o results
fastnlr -i genome.fasta -o results -t 8
fastnlr -i genome.fasta -o results --checkpoint ckpt/ --flexible-seed 3
fastnlr -i genome.fasta --gff annotation.gff3 -o annotation_results
```

`-i` and `-o` are required; the built-in motif profiles are used unless `-x` / `-y` are given.
With `--gff`, `-i` remains a genome FASTA and the annotation is used to translate CDS.
Text protein FASTA input is not supported.
Direct protein FASTA input is not supported.

### Flag groups

| Group | Flags |
|---|---|
| Input | `-i <INPUT>` |
| Motif config (genome mode only) | `-x <mot.txt>`, `-y <store.txt>` |
| Output | `-o <DIR>`, `--flank <bp>` (default 2000) |
| Performance | `-t <threads>`, `-n <fragments per chunk>` (default 1000) |
| Run control | `--checkpoint <DIR>` |
| Observability | `--progress <auto\|bar\|simple\|off>`, `--log-level <level>` |
| Recall (genome mode only) | `--motif-accept-p`, `--motif-prelim-p`, `--relaxed-seed <n>`, `--flexible-seed <n>`, `--require-ploop` |
| Motif library metadata | `--motif-category <ID=CAT>`, `--seed-combination <IDS>`, `--signature <IDS>` |
| GFF3 mode | `--gff <GFF3>`, `--hmm <FILE>` (repeatable) |

## GFF3 mode

`--gff` replaces six-frame motif scanning with GFF3-guided CDS translation:

- **Transcript selection** — genes are parsed from `gene → mRNA/transcript → CDS`; the transcript
  with the largest total CDS length is selected for each gene. Gzipped GFF3 is accepted.
- **Translation** — CDS segments are concatenated in transcript order and reverse-complemented on
  the negative strand. The first CDS phase is applied as a 5' translation offset while the raw
  CDS FASTA is preserved. Trailing partial codons are trimmed and recorded. A terminal stop is
  removed; internal stops become `X` and are counted.
- **Domains** — 23 Pfam models are embedded and scored with their gathering thresholds. Five
  models identify NB-ARC, TIR and N-terminal CC families; 18 models cover generic, plant-specific
  and terminal LRR regions. LRR scanning runs only after NB-ARC or TIR passes GA.
- **Classification** — `nterm` is the N-terminal-most TIR/CC domain before the NB-ARC domain;
  `class` collapses adjacent equal categories using their real protein coordinates. LRR before
  NB-ARC remains diagnostic evidence but is not C-terminal LRR;
  `has_lrr_any` reports all passed LRR, `has_cterm_lrr` requires LRR after NB-ARC, and
  `complete` means NB-ARC + C-terminal LRR;
  `family` distinguishes TNL / RNL / CNL / NB-only / TIR-only (RPW8 CC is RNL, non-RPW8 CC is
  CNL); a protein is an NLR candidate when it has NB-ARC, or TIR together with C-terminal LRR.
- **PWM and DNA-only flags are rejected** in this mode: `-x`, `-y`, `--motif-accept-p`,
  `--motif-prelim-p`, `--flank`, `--relaxed-seed`,
  `--flexible-seed`, `--require-ploop`, `--signature`, `--seed-combination`, `--checkpoint`,
  `--tir-accept-p` and partial-locus filtering.

Outputs: `<prefix>.gff.cds.fasta`, `.gff.proteins.fasta`, `.gff.domains.tsv`,
`.gff.genes.tsv`, `.gff.nlr.tsv`, `.gff.nlr.gff3`, `.gff.skipped.tsv`,
`.gff.summary.txt` and `.gff.summary.json`. The six-frame locus outputs are not produced.

## Pipeline

`run` / `run_with_progress` execute the following stages:

1. **Load the profile and build the rule definition.** User paths take precedence over the
   embedded defaults; CLI-declared categories, seeds, signatures and the flexible-seeding
   configuration are injected into `AnnotatorSignatureDefinition`.
2. **Strict checkpoint short-circuit.** `manifest.json` fingerprints the input, profiles,
   thresholds, slicing/assembly parameters and external rules. Missing, incomplete, corrupt or
   mismatched checkpoints fail explicitly. The first-run timestamp is retained so direct and
   resumed runs produce byte-identical `.nlr.txt`, `.nlr.gff` and `.nlr.bed` files.
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
6. **Write the checkpoint atomically** after a complete, non-interrupted scan.

`SIGINT` sets an atomic flag; the scan stops after the current chunk and the completed results
are still written out.

## Public API

| Item | Description |
|---|---|
| `RunConfig` | full run configuration, see below |
| `RunResult` | `nlrs`, `motifs_by_seq`, `def`, and an optional checkpoint-fixed GFF timestamp |
| `run(&RunConfig)` / `run_with_progress(&RunConfig, Option<&ProgressBar>)` | pipeline entry points; the second reports progress |
| `all_motifs(&RunResult)` | flattens motif hits across contigs |
| checkpoint internals | versioned `manifest.json` plus strict motif TSV import |
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
| `flexible_seed` | `None` in the library constructor; CLI default is `min_nbarc = 3` | `--flexible-seed` / `--require-ploop`; `--no-flexible-seed` disables it |

## Outputs

For an input `genome.fasta`, the prefix `genome` is derived by stripping `.fa/.fasta/.fna/.fas`
and `.gz`, and the following files are written into `-o`:

| File | Writer |
|---|---|
| `genome.nlr.txt` | `nlr-output` |
| `genome.nlr.gff` (GFF2) | `nlr-output` |
| `genome.nlr.bed` | `nlr-output` |
| `genome.motifs.bed` | `nlr-output` |
| `genome.nbarc.fasta` | `nlr-output` |
| `genome.loci.fasta` | `nlr-output` |
| `genome.stats.tsv` | `nlr-report` |
| `genome.summary.txt` | this crate |

All outputs are written by default; there is no opt-in flag per format.
PNG statistics plots are not generated.

GFF3 mode instead writes the `.gff.*` file set described in [GFF3 mode](#gff3-mode).

## Notes

- Checkpoints are reused only with identical input/profile hashes and runtime settings; the saved
  motif set is assembled with the current validated `AssembleParams`.
- Locus extraction streams the genome contig by contig: each contig is read once, used to append
  its locus FASTA records and to annotate the HMM N-terminal domain, then dropped. Peak memory in
  the output stage no longer scales with the whole genome.
- The GFF date header uses the host's local timezone; checkpoint runs retain the first-run
  timestamp so resumed output is reproducible.
- Peak memory is reported on Linux (`/proc/self/status`) and macOS (`getrusage`); on other
  platforms the summary prints `n/a` instead of a misleading `0.00 MB`.

## Testing

```bash
cargo test -p fastnlr
```

`tests/run.rs` drives the library pipeline from synthetic DNA and asserts that motifs are
recovered; `tests/e2e.rs` covers the scan-to-assembly path on back-translated constructs.
`tests/outputs.rs` verifies the main output files and BED block structure;
`tests/annotation_mode.rs` covers GFF3 translation, HMM-only classification, deterministic
thread/batch output, pre-NBARC LRR rejection and rejection of PWM/genome-only flags.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
