<div align="center">

# FastNLR

**Fast, self-contained NLR immune-receptor locus annotation for plant genomes**

A high-speed, accurate NLR annotation tool written in Rust for plant genomes

[![bioconda](https://anaconda.org/bioconda/fastnlr/badges/version.svg)](https://anaconda.org/bioconda/fastnlr)
[![Downloads](https://anaconda.org/bioconda/fastnlr/badges/downloads.svg)](https://anaconda.org/bioconda/fastnlr)
[![Platforms](https://anaconda.org/bioconda/fastnlr/badges/platforms.svg)](https://anaconda.org/bioconda/fastnlr)
[![CI](https://github.com/CropCoder/FastNLR/actions/workflows/rust.yml/badge.svg)](https://github.com/CropCoder/FastNLR/actions)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-1.5.2-green.svg)](https://github.com/CropCoder/FastNLR/releases)

[Install](#install) · [Quick start](#quick-start) · [Features](#features) · [Performance](#performance) · [Documentation](#documentation) · [Citation](#citation)

</div>

---

FastNLR scans the six reading frames of a genome assembly for amino-acid **motifs** and
assembles them into **NLR** (Nucleotide-binding, Leucine-rich Repeat) immune-receptor loci.
It combines multithreading, memory-mapped FASTA input, embedded default motif profiles,
checkpoint resume, and coordinate-consistent multi-format outputs.

NLR genes encode a major class of plant intracellular immune receptors and are a key target of
disease-resistance breeding. FastNLR identifies complete and partial NLR loci directly from an
assembled genome, or classifies NLR domains from a protein FASTA. It is distributed through
[bioconda](https://anaconda.org/bioconda/fastnlr) as well as prebuilt release binaries and source.

## Install

### Recommended — bioconda

```bash
conda install bioconda::fastnlr
```

The `bioconda::` prefix pins the package to the bioconda channel. If the channels are not already
in your conda configuration, give them explicitly (`conda-forge` must come first):

```bash
conda install -c conda-forge -c bioconda bioconda::fastnlr
# mamba / micromamba are drop-in and much faster
mamba install -c conda-forge -c bioconda bioconda::fastnlr
```

### Prebuilt binary

```bash
curl -L -o fastnlr.tar.gz \
  https://github.com/CropCoder/FastNLR/releases/download/v1.5.2/fastnlr-x86_64-unknown-linux-gnu.tar.gz
tar -xzf fastnlr.tar.gz
./fastnlr --version
```

Release assets follow the `fastnlr-<target>.tar.gz` naming scheme (Windows builds use `.zip`,
macOS assets are also `.tar.gz`) and ship with a SHA-256 checksum file.

### From source

```bash
# Requires Rust 1.70+ (developed on 1.94.1)
git clone https://github.com/CropCoder/FastNLR.git
cd FastNLR
cargo build --release
# binary: target/release/fastnlr
```

Package page: <https://anaconda.org/bioconda/fastnlr> · Packaging recipe: [`conda/`](conda/)

## Quick start

```bash
# 1. Genome (DNA) mode — built-in motif library, 8 threads
#    -> results/<genome>.nlr.txt, .nlr.gff, .nlr.bed, .motifs.bed, .nbarc.fasta,
#       .loci.fasta, .stats.tsv, .summary.txt, .plots/
fastnlr -i genome.fasta -o results -t 8

# 2. Protein mode — HMMER domain scan on a protein FASTA (no genome coordinates)
#    -> protein_results/<prefix>.domains.tsv, .nlr.tsv, .summary.txt, .plots/
fastnlr --protein -i proteins.fa -o protein_results

# 3. Resume a run (skips the expensive scan)
fastnlr -i genome.fasta -o results --checkpoint ckpt/
```

Run `fastnlr --help` for the full flag reference, or see the
[project page](https://cropcoder.github.io/FastNLR/) for annotated examples, output-format
details and benchmark results.

## Features

- **Zero-config, self-contained binary** — the standard `mot.txt` (PWM) and `store.txt` (CDF)
  configs are embedded at compile time; override with `-x`/`-y` for custom motifs.
- **High performance** — Rust + rayon multithreading, memory-mapped FASTA parsing and SIMD
  scoring; throughput holds at roughly 3.7 Mb/s across two orders of magnitude of genome size.
- **Flexible seeding, on by default** — category/rank-based seeding runs out of the box and
  recovers divergent NLR/RNL loci. It is **additive**: it never drops a locus that the exact
  motif-id combinations would have found (`--no-flexible-seed` restores the strict behaviour).
- **HMM-confirmed domain architecture** — each assembled locus is re-scanned with the same
  bundled HMMER profiles used in protein mode, so N-terminal CC/TIR domains that the PWM library
  misses are still reported instead of collapsing to a bare `NBARC-LRR`.
- **Tandem-repeat aware assembly** — repeated or duplicated domain units and trailing NB-ARC
  segments are merged back into their parent locus, so one gene is not reported as several loci.
- **Extensible motif libraries** — an optional RNL/helper-CC + TIR library (28 motifs) recovers
  RNL/helper NLR and TIR-only loci that the built-in 20-motif profile misses.
- **Protein input mode** — `--protein` classifies domain architecture with a bundled HMMER 3.4
  port (Pfam NB-ARC / TIR / Rx_N / RPW8 / R1-type at gathering thresholds) plus the PWM LRR motifs.
- **Resumable runs, rich reporting, graceful interruption** — `--checkpoint` skips the scan on
  reruns; every run writes a human-readable summary, per-chromosome and per-motif TSV statistics
  and PNG plots; Ctrl-C still emits the batches already completed.

## Performance

Benchmarked on 18 plant genomes from 197 Mb to 4.7 Gb, 8 threads, warm cache, median of 3 runs:

| Genome | Size | Wall time | Peak RSS |
|--------|------|-----------|----------|
| *Cucumis sativus* (cucumber) | 197 Mb | 53.5 s | 210 MB |
| *Glycine max* (soybean) | 991 Mb | 269.5 s | 968 MB |
| *Hordeum vulgare* (barley) | 4.3 Gb | 1174 s | 4159 MB |
| *Zanthoxylum bungeanum* (pepper) | 4.7 Gb | 1230 s | 4625 MB |

Peak memory stays close to 1:1 with genome size (ratios 0.97–1.15), so budget RAM of the same
order as the assembly. Multithreading scales 13–19× from 4 to 128 threads on the largest genomes;
keep `-t` at or below the physical core count. The `sensitive` recall preset
(`--motif-accept-p`/`--motif-prelim-p` relaxed) is diagnostic only: it inflates runtime by
roughly 40× and is not intended for production runs.

**Accuracy** (RefPlantNLR, 400 curated NLR CDS): **99.75%** gene recall and **96.49%**
domain-architecture agreement — rising to **98.19%** when the 12 reference entries whose
`Structure` and `Subclass` columns contradict each other are excluded. Per-domain F1 is
1.000 (NB-ARC), 0.999 (LRR), 0.995 (TIR) and 0.989 (CC).

## Documentation

| Topic | Where |
|---|---|
| Full flag reference, output formats, worked examples | <https://cropcoder.github.io/FastNLR/> |
| Extended RNL/TIR motif library and design notes | [`data/motif_library_rnl_tir/`](data/motif_library_rnl_tir/) |
| Flexible-seeding upgrade design | [`docs/flexible-seeding-upgrade.md`](docs/flexible-seeding-upgrade.md) |
| Packaging recipe | [`conda/`](conda/) |
| Development log, benchmarks and design records | [`Runlog.md`](Runlog.md) |

## Development

```bash
cargo build --release                              # optimized build
cargo test --workspace                             # 69 tests
cargo fmt --all -- --check                         # formatting
cargo clippy --workspace --all-targets -- -D warnings
```

The release profile uses `opt-level=3`, `lto="fat"`, `codegen-units=1` and `panic="abort"`.
CI (`.github/workflows/rust.yml`) builds and tests on every push; pushing a `vX.Y.Z` tag builds
Linux/Windows/macOS binaries and uploads them to GitHub Releases.

## Citation

FastNLR is developed and maintained by Jiwen Zhao. A dedicated manuscript is in preparation;
until then, please cite this repository:

> Jiwen Zhao. FastNLR: a high-speed, accurate NLR annotation tool. https://github.com/CropCoder/FastNLR

## License

GPL-3.0-only. See [LICENSE](LICENSE).

## Author

**Jiwen Zhao** — <https://github.com/CropCoder> · Contributions and issue reports are welcome.
