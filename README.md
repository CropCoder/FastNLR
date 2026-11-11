<div align="center">

<img src="data/motif_library_rnl_tir/logo.png" alt="FastNLR logo" width="640">

**Fast, self-contained NLR immune-receptor locus annotation for plant genomes**

[![bioconda](https://anaconda.org/bioconda/fastnlr/badges/version.svg)](https://anaconda.org/bioconda/fastnlr)
[![Downloads](https://anaconda.org/bioconda/fastnlr/badges/downloads.svg)](https://anaconda.org/bioconda/fastnlr)
[![CI](https://github.com/CropCoder/FastNLR/actions/workflows/rust.yml/badge.svg)](https://github.com/CropCoder/FastNLR/actions)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.85%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-1.9.1-green.svg)](https://github.com/CropCoder/FastNLR/releases)

[Install](#install) · [Quick start](#quick-start) · [Features](#features) · [Performance](#performance) ·
[Interfaces](#other-interfaces) · [Website](https://cropcoder.github.io/FastNLR/) · [Citation](#citation)

</div>

---

FastNLR scans the six reading frames of a genome assembly for amino-acid **motifs** and
assembles them into **NLR** (Nucleotide-binding, Leucine-rich Repeat) immune-receptor loci.
It combines multithreading, memory-mapped FASTA input, embedded default motif profiles,
checkpoint resume, and coordinate-consistent multi-format outputs.

NLR genes encode a major class of plant intracellular immune receptors and are a key target of
disease-resistance breeding. FastNLR identifies complete and partial NLR loci directly from an
assembled genome, or uses an optional GFF3 annotation to translate and classify gene proteins. It is distributed through
[bioconda](https://anaconda.org/bioconda/fastnlr) as well as prebuilt release binaries and source.

## Install

### Recommended — bioconda

```bash
conda install bioconda::fastnlr
```

### Prebuilt binary

```bash
curl -L -o fastnlr-v1.9.1-x86_64-unknown-linux-gnu.tar.gz \
  https://github.com/CropCoder/FastNLR/releases/download/v1.9.1/fastnlr-v1.9.1-x86_64-unknown-linux-gnu.tar.gz
tar -xzf fastnlr-v1.9.1-x86_64-unknown-linux-gnu.tar.gz
./fastnlr --version
```

Release assets follow the `fastnlr-<tag>-<target>.tar.gz` naming scheme (Windows builds use
`.zip`, macOS assets are also `.tar.gz`). Linux x86_64 binaries target glibc 2.17, so they run
on CentOS 7, RHEL 8, Debian 11, Ubuntu 20.04 and newer; on older systems use the bioconda
package or build from source.

### Desktop application — Windows and macOS

FastNLR desktop builds provide genome-only and GFF3-guided workflows in a native application. Download
them from the latest [GitHub Release](https://github.com/CropCoder/FastNLR/releases/latest).

| Platform | Download | Requirements |
|---|---|---|
| Windows GUI | NSIS `.exe` installer | Windows 10/11 x64, WebView2 |
| macOS GUI | `.dmg` installer | macOS 12+, Apple Silicon |
| Linux CLI | x86_64 `.tar.gz` | glibc 2.17+ |

The initial desktop builds are not yet signed with a commercial Windows certificate or notarized
with an Apple Developer ID. Windows SmartScreen and macOS Gatekeeper may therefore show a warning;
these packages are suitable for evaluation while signing credentials are being prepared.

### From source

```bash
# Requires Rust 1.85+
git clone https://github.com/CropCoder/FastNLR.git
cd FastNLR
cargo build --release --locked
# binary: target/release/fastnlr
```

Package page: <https://anaconda.org/bioconda/fastnlr> · Packaging recipe: [`conda/`](conda/)

## Quick start

```bash
# 1. Genome (DNA) mode — built-in motif library, 8 threads
#    -> results/<genome>.nlr.txt, .nlr.gff, .nlr.bed, .motifs.bed, .nbarc.fasta,
#       .loci.fasta, .stats.tsv, .summary.txt, .summary.json
fastnlr -i genome.fasta -o results -t 8

# 2. GFF3 mode — translate the longest CDS transcript per annotated gene,
#    then classify the proteins with ordered HMM evidence only.
#    -> annotation_results/<genome>.gff.{cds.fasta,proteins.fasta,domains.tsv,
#       genes.tsv,nlr.tsv,nlr.gff3,skipped.tsv,summary.txt,summary.json}
fastnlr -i genome.fasta --gff annotation.gff3 -o annotation_results

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
  scoring; throughput holds at roughly 3.7 Mb/s across the benchmarked 197 Mb–4.7 Gb range.
- **Flexible seeding, on by default** — category/rank-based seeding runs out of the box and
  recovers divergent NLR/RNL loci. It is **additive**: it never drops a locus that the exact
  motif-id combinations would have found (`--no-flexible-seed` restores the strict behaviour).
- **HMM-assisted N-terminal domains** — when the PWM library does not identify a CC/TIR domain,
  the locus is checked with the same bundled HMMER profiles used in GFF3-guided protein mode. HMM hits fill
  that missing N-terminal domain instead of leaving a bare `NBARC-LRR`.
- **Tandem-repeat aware assembly** — repeated or duplicated domain units and trailing NB-ARC
  segments are merged back into their parent locus, so one gene is not reported as several loci.
- **High-precision defaults** — TIR motifs use a stricter acceptance threshold
  (`--tir-accept-p`, default `1e-7`), and post-assembly filtering removes unsupported partial
  loci. Use `--keep-partial-loci` to restore legacy behaviour.
- **RNL-aware default library** — the default profile includes the classic 20 motifs plus
  RNL/helper-CC profiles 21 and 23 and conservative `CC,1` seeds. This recovers ADR1-L1
  without relaxing the global PWM p-value threshold. An optional 28-motif RNL/TIR library
  remains available for broader exploratory scans.
- **GFF3-guided protein mode** — with `--gff`, FastNLR selects the longest CDS transcript per
  annotated gene, streams the genome to extract and translate CDS, then classifies the proteins
  with 23 bundled Pfam HMMs. The 18 LRR models are scanned only for NB-ARC/TIR candidates.
  LRR must lie after NB-ARC to count as C-terminal evidence. The output reports
  gene/transcript coordinates, protein/CDS/genomic domain coordinates, structural class,
  completeness and biological family (`TNL` / `RNL` / `CNL` / `NB-only` / `TIR-only`).
- **Genome-only remainder** — direct protein FASTA input is no longer supported. Without GFF3,
  FastNLR runs the original six-frame motif assembly pipeline.
- **Strict resumable runs** — checkpoints include input/profile/configuration fingerprints and are
  reused only when the current run is compatible. Incomplete or stale checkpoints fail explicitly;
  every run writes human-readable and JSON summaries plus per-chromosome/per-motif TSV statistics.

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

**TAIR12 reference validation:** FastNLR 1.7.2 completed the standard DNA-mode workflow on
*Arabidopsis thaliana* TAIR12 (`GCA_978657495.1`, 142.5 Mb) with the built-in motif profiles
and 8 threads. It scanned 7,918 fragments, assembled 4,929 motif hits into 210 NLR loci, and
finished in 43.6 s with 162 MB peak RSS.

**Accuracy** (RefPlantNLR, 400 curated NLR CDS): **99.75%** gene recall and **96.49%**
domain-architecture agreement — rising to **98.19%** when the 12 reference entries whose
`Structure` and `Subclass` columns contradict each other are excluded. Per-domain F1 is
1.000 (NB-ARC), 0.999 (LRR), 0.995 (TIR) and 0.989 (CC).

## Other interfaces

- **Web server** — `fastnlr-server` is a self-contained axum service with an embedded UI and a
  job-oriented REST API. Build with `cargo build --release -p nlr-server`, then run
  `fastnlr-server --port 8787`. It is also packaged as a local TBtools web plugin via
  `packaging/Build-TBtools-Plugin-V1.sh`.
- **Desktop app** — the Tauri 2 application in `app/` provides native Windows and macOS workflows,
  live progress and cancellation, persistent run history, result tables, summary charts, and an
  interactive chromosome panorama with locus filtering and density curves while reusing the same
  `nlr-service` execution layer.
- **AI agents** — the `skills/fastnlr/` directory provides a portable Agent Skill for Codex and
  Claude Code.
- **Publication plots** — the `visualization/` R scripts turn `summary.json` into chromosome maps
  and statistics PDFs plus an HTML report.

## Development

```bash
cargo build --release --locked                     # optimized build
cargo test --workspace --locked                    # run all Rust tests
cargo fmt --all -- --check                         # formatting
cargo clippy --workspace --all-targets -- -D warnings
```

Desktop development requires Node.js 22+ and Rust 1.90+:

```bash
cd app
npm ci
npm run typecheck
npm test
npm run tauri build
```

The release profile uses `opt-level=3`, `lto="fat"`, `codegen-units=1` and `panic="abort"`.
CI (`.github/workflows/rust.yml`) builds and tests pushes to `main` and pull requests; pushing a
`vX.Y.Z` tag builds the Linux CLI plus Windows/macOS GUI installers and uploads them to GitHub
Releases.

## Citation

A dedicated manuscript is in preparation. Until then, please cite this repository:

> FastNLR: a high-speed, accurate NLR annotation tool. https://github.com/CropCoder/FastNLR

## License

GPL-3.0-only. See [LICENSE](LICENSE).
