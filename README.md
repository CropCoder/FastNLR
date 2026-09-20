<div align="center">

# FastNLR

**Fast, self-contained NLR immune-receptor locus annotation for plant genomes**

A high-speed, accurate NLR annotation tool written in Rust for plant genomes

[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![CI](https://github.com/CropCoder/FastNLR/actions/workflows/rust.yml/badge.svg)](https://github.com/CropCoder/FastNLR/actions)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-1.1.1-green.svg)](https://github.com/CropCoder/FastNLR/releases)
[![Platform](https://img.shields.io/badge/platform-linux%20x86__64-lightgrey.svg)](https://github.com/CropCoder/FastNLR/releases)

[Features](#features) · [Quick Start](#quick-start) · [Usage](#usage) · [Extended motif library](#extended-motif-library) · [Output Formats](#output-formats) · [Architecture](#architecture) · [Citation](#citation)

</div>

---

FastNLR scans the six reading frames of a genome assembly for amino-acid **motifs** and assembles them into **NLR** (Nucleotide-binding, Leucine-rich Repeat) immune-receptor loci. It combines multithreading, memory-mapped FASTA input, embedded default motif profiles, checkpoint resume, and coordinate-consistent multi-format outputs.

NLR genes encode a major class of plant intracellular immune receptors and are a key target of disease-resistance breeding. FastNLR provides a self-contained pipeline for identifying complete and partial NLR loci directly from an assembled genome. The project is under active development, with future releases targeting publication-grade documentation, validation, and packaging.

## Features

- **Zero-config, self-contained binary** — the standard `mot.txt` (PWM) and `store.txt` (CDF) configs are embedded at compile time. Run on any FASTA with no extra files; override with `-x`/`-y` when you need custom motifs.
- **Extensible motif libraries** — supports more than the built-in 20 motifs, plus CLI-declared domain categories, seed combinations, and signatures. An optional RNL/helper-CC + TIR library is included for recovering RNL and TIR-only loci.
- **High performance** — Rust + [rayon](https://github.com/rayon-rs/rayon) multithreading, memory-mapped large-file FASTA parsing, and SIMD (`wide`) cross-window scoring. Typical plant genomes finish in seconds to minutes.
- **Coordinate-consistent output** — `-f` loci extraction matches the GFF/BED coordinates exactly; see [Implementation notes](#implementation-notes).
- **Resumable runs** — `--checkpoint` saves motif results after the scan; reruns skip the expensive scan and go straight to assembly.
- **Rich reporting** — human-readable summary, TSV statistics (global / per-chromosome / per-motif), and PNG plots out of the box.
- **Graceful interruption** — Ctrl-C outputs whatever batches have completed instead of dropping everything.
- **Familiar flag interface** — supports `-i/-x/-y/-o/-g/-b/-m/-a/-f/-c/-t/-n` workflows.

## Quick Start

### Option A — download the prebuilt binary

```bash
# Download the Linux x86_64 asset from Releases, then:
tar -xzf fastnlr-x86_64-unknown-linux-gnu.tar.gz
./fastnlr --version
./fastnlr -i genome.fasta -o result.txt -p result
```

Verify integrity:

```bash
sha256sum -c fastnlr-x86_64-unknown-linux-gnu.tar.gz.sha256
```

Release assets use the `fastnlr-<target>.tar.gz` naming scheme; Windows builds use `.zip`, and macOS assets are also `.tar.gz`.

### Option B — build from source

```bash
# Requires Rust 1.70+ (developed on 1.94.1)
git clone https://github.com/CropCoder/FastNLR.git
cd FastNLR
cargo build --release
# binary: target/release/fastnlr
```

## Usage

```bash
fastnlr -i <input.fasta> -o <output.txt> [output flags] [options]
```

### Examples

```bash
# 1. Basic loci annotation (uses built-in mot.txt / store.txt)
fastnlr -i genome.fasta -o out.txt -g out.gff -b out.bed

# 2. Prefix-derived subfiles + multithreading
#    -> out.nlr.txt, out.nlr.gff, out.nlr.bed, out.motifs.bed, out.nbarc.fasta
fastnlr -i genome.fasta -o out.txt -p out -t 8

# 3. Full run: report + summary + plots + checkpoint resume
fastnlr -i genome.fasta -o out.txt -p out \
  --stats stats.tsv --summary --plot plots/ --checkpoint ckpt/

# 4. Extract NLR loci sequences (±2000 bp flanking) across ALL contigs
fastnlr -i genome.fasta -o out.txt -p out -f genome.fasta loci.fasta 2000

# 5. Custom motif config override
fastnlr -i genome.fasta -x custom_mot.txt -y custom_store.txt -o out.txt -p out

# 6. Resume from checkpoint after an interrupted run
fastnlr -i genome.fasta -o out.txt -p out --checkpoint ckpt/

# 7. Optional RNL/helper-CC + TIR extended library (28 motifs)
source data/motif_library_rnl_tir/flags.sh
fastnlr -i genome.fasta -o out.txt -p out \
  -x data/motif_library_rnl_tir/mot.txt \
  -y data/motif_library_rnl_tir/store.txt \
  --motif-category "$TIS_MOTIF_CATEGORY" \
  --seed-combination "$TIS_SEED_COMBINATION" \
  --signature "$TIS_SIGNATURE"
```

### Flags

**Input**

| Flag | Description |
|------|-------------|
| `-i <fasta>` | Input genome FASTA (may be gzip-compressed). **Required.** |
| `-x <mot.txt>` | PWM config (optional; default: built-in). |
| `-y <store.txt>` | CDF config (optional; default: built-in). |

**Motif library metadata** — only needed for custom libraries with more than 20 motifs.

| Flag | Description |
|------|-------------|
| `--motif-category ID=CAT` | Declare the domain category for a motif id outside the built-in tables, e.g. `21=CC,25=TIR`. Repeatable or comma-separated. |
| `--seed-combination IDS` | Add a seed motif-id combination, e.g. `21,4`; multiple combinations are separated by `;`. |
| `--signature IDS` | Add a fragment signature, e.g. `25,26`; multiple signatures are separated by `;`. |

**Sensitivity / recall** — optional tuning for divergent NLRs.

| Flag | Description |
|------|-------------|
| `--motif-accept-p <p>` | Final motif-accept p-value threshold (default `1e-5`). |
| `--motif-prelim-p <p>` | Fragment prefilter p-value threshold (default `1e-4`). |
| `--relaxed-seed <n>` | Relaxed seeding: accept a consecutive hit run with at least `n` NB-ARC motifs (off by default). |

**Output** — `-o` is required; add any other format as needed.

| Flag | Description |
|------|-------------|
| `-o <txt>` | NLR loci report (tabular). **Required.** |
| `-g <gff>` | NLR loci (GFF3). |
| `-b <bed>` | NLR loci (BED12, color-coded). |
| `-m <bed>` | Motif intervals (BED). |
| `-a <fasta>` | NB-ARC multiple alignment (fasta). |
| `-f <genome> <out> <bp>` | Loci sequences (fasta) with flanking bp; extracts from every contig. |
| `-c <tsv>` | Export precomputed motif results (reusable as checkpoint import). |
| `-p, --output-prefix <p>` | Auto-derive all of the above as `p.nlr.txt`, `p.nlr.gff`, … |

**Run control**

| Flag | Description |
|------|-------------|
| `-t <n>` | Threads (default: auto-detect). |
| `-n <n>` | Fragments per batch (default 1000). |
| `--checkpoint <dir>` | Save/load motif results to skip rescan. |
| `--tmpdir <dir>` | Temp directory. |
| `--progress <auto\|bar\|simple\|off>` | Progress bar (default auto). |
| `--stats <file>` | TSV statistics report. |
| `--plot <dir>` | PNG statistics plots. |
| `--summary` | Per-chromosome summary to stdout. |
| `--log-level <lvl>` | trace/debug/info/warn/error (default info). |

Run `fastnlr --help` for the full, grouped reference.

## Extended motif library

The default build ships the built-in 20-motif profile. The repository also contains an optional library at [`data/motif_library_rnl_tir/`](data/motif_library_rnl_tir/) that adds eight motifs:

- `21`–`24`: RNL / helper-NLR N-terminal CC motifs (`ADR1`, `NRG1`, and Solanaceous `NRC`).
- `25`–`28`: TIR-family motifs, including short and TIR-only loci that the standard `[18,15,13]` seed could miss.

Because these motifs live outside the built-in rule tables, enable them together with their declared categories, seeds, and signatures:

```bash
source data/motif_library_rnl_tir/flags.sh
fastnlr -i genome.fasta -o out.txt -p out \
  -x data/motif_library_rnl_tir/mot.txt \
  -y data/motif_library_rnl_tir/store.txt \
  --motif-category "$TIS_MOTIF_CATEGORY" \
  --seed-combination "$TIS_SEED_COMBINATION" \
  --signature "$TIS_SIGNATURE"
```

`flags.sh` contains the complete parameter lists; adjust `-x`/`-y` paths if you copy the files elsewhere. The companion [`data/motif_library_rnl_tir/README.md`](data/motif_library_rnl_tir/README.md) records the motif-design rationale, benchmark results, and regeneration workflow.

For sensitivity tuning on divergent NLRs, three additional flags are available: `--motif-accept-p`, `--motif-prelim-p`, and `--relaxed-seed` (see the flag tables below).

## Output Formats

- **`-o` report** — one row per NLR: `seqname, name, domain-class, start, end, strand, motif-list`.
- **`-g` GFF3** — header relabeled for FastNLR with a live system timestamp; `source` column = `FastNLR`.
  ```
  ##gff-version 2
  ##source-version FastNLR V1.1.1
  ##date 2026-08-20 09:30:01
  ##Type DNA
  ```
- **`-b` BED12** — color-coded blocks (green = complete, orange = partial, red = contains stop codon); reverse-strand blocks reversed.
- **`-a` NB-ARC alignment** — P-loop-anchored multiple alignment with gap-padding for missing motifs.
- **`-f` loci fasta** — extracted ±flanking sequence, reverse-complemented on the reverse strand, 100 bp per line.

## Architecture

FastNLR is a layered Cargo workspace — each crate has a single responsibility and the dependency graph flows bottom-up:

```
nlr-core     domain model + rule tables (rank/class/color/seed/signature/consensus; built-in + CLI-declared extras)
nlr-config   parse mot.txt (PWM) / store.txt (CDF); built-in embed and >20-motif library support
nlr-seq      six-frame translation, reverse complement, codon table, FASTA reader, chopper
nlr-scan     sliding-window scoring + non-overlap arbitration + signature pre-filter (SIMD)
nlr-assemble findSeeds -> mergeSeeds -> elongate three-step assembly
nlr-output   txt/GFF/BED/motifBED/alignment-fasta/loci-fasta/TSV output
nlr-report   run statistics aggregation
nlr-plot     plotters statistics plots (bundled font backend, no system font deps)
nlr-cli      clap CLI entry + pipeline orchestration (rayon + checkpoint + SIGINT); binary: fastnlr
```

**Pipeline:** `FASTA → chop (overlap) → six-frame translate → scan (PWM+CDF) → signature filter → coordinate map → three-step assembly → multi-format output`.

## Implementation notes

FastNLR enforces coordinate-consistent behavior across its output formats:

1. `-a` uses corrected P-loop location logic.
2. `-a` replaces stop codons and unknown amino acids with `_`.
3. `-f` applies a unified coordinate clamp on the last contig.
4. `-f` extracts loci from every contig, not only the first.
5. `-f` keeps extracted sequence coordinates consistent with the reported GFF/BED coordinates.
6. `##date` reflects the current system time.
7. GFF `##source-version` and the `source` column are labeled `FastNLR`.
8. p-values in motif BED and export TSV use scientific notation for very small values.

## Correctness

The pipeline is validated by unit and integration tests covering coordinate mapping (both strands), six-frame translation tail offsets, the codon table, seed/signature tables, motif scanning, three-step assembly, and every output format.

```bash
cargo test --workspace      # 40 tests, all passing
```

## Performance

| Genome | Size | Threads | Wall time ( indicative ) |
|--------|------|---------|--------------------------|
| *Arabidopsis* chr1 | ~30 Mb | 8 | a few seconds |
| Rice genome | ~370 Mb | 16 | tens of seconds |
| Wheat genome | ~14 Gb | 32 | minutes |

Memory: memory-mapped FASTA keeps peak RAM roughly proportional to `threads × batch_size × fragment_length`, independently of genome size.

## Development

```bash
cargo build                          # debug build
cargo build --release                # optimized build
cargo test --workspace               # run all tests
cargo bench                          # benchmarks (if enabled)
```

The release profile uses `opt-level=3`, `lto="fat"`, `codegen-units=1`, and `panic="abort"` for maximum performance and a small binary.

Continuous integration and release automation are defined in `.github/workflows/`:

- `rust.yml` runs workspace build and tests on pushes and pull requests to `main`.
- `release.yml` builds `fastnlr` for Linux, Windows, and macOS x86_64 when a `vX.Y.Z` tag is pushed, then uploads archives and SHA-256 checksums to GitHub Releases.

The static project page is generated from [`docs/index.html`](docs/index.html) and published through GitHub Pages.

## Citation

FastNLR is developed and maintained by Jiwen Zhao. A dedicated manuscript is in preparation; until then, please cite this repository:

> Jiwen Zhao. FastNLR: a high-speed, accurate NLR annotation tool. https://github.com/CropCoder/FastNLR

## License

GPL-3.0-only. See [LICENSE](LICENSE).

## Author

**Jiwen Zhao** — https://github.com/CropCoder

Contributions and issue reports are welcome.
