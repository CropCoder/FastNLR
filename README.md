<div align="center">

# FastNLR

**Fast, self-contained NLR immune-receptor locus annotation for plant genomes**

A high-performance Rust rewrite of [NLR-Annotator](https://pubmed.ncbi.nlm.nih.gov/32184345/)

[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org)
[![Version](https://img.shields.io/badge/version-1.0.0-green.svg)](https://github.com/CropCoder/FastNLR/releases)
[![Platform](https://img.shields.io/badge/platform-linux%20x86__64-lightgrey.svg)](https://github.com/CropCoder/FastNLR/releases)

[Features](#features) · [Quick Start](#quick-start) · [Usage](#usage) · [Output Formats](#output-formats) · [Architecture](#architecture) · [Citation](#citation)

</div>

---

FastNLR scans the six reading frames of a genome assembly for amino-acid **motifs** and assembles them into **NLR** (Nucleotide-binding, Leucine-rich Repeat) immune-receptor loci. It is a from-scratch Rust port of the Java tool *NLR-Annotator* (Steuernagel et al., 2020), preserving the original algorithm semantics while adding multithreading, embedded default configs, checkpoint resume, and several bug fixes.

NLR genes encode the largest class of plant intracellular immune receptors and are a key target of disease-resistance breeding. FastNLR lets you survey any assembled genome for complete and partial NLR loci in minutes.

## Features

- **Zero-config, self-contained binary** — the standard `mot.txt` (PWM) and `store.txt` (CDF) configs are embedded at compile time. Run on any FASTA with no extra files; override with `-x`/`-y` when you need custom motifs.
- **High performance** — Rust + [rayon](https://github.com/rayon-rs/rayon) multithreading, memory-mapped large-file FASTA parsing, and SIMD (`wide`) cross-window scoring. Typical plant genomes finish in seconds to minutes.
- **Coordinate-consistent output** — `-f` loci extraction matches the GFF/BED coordinates exactly. (The original Java `writeNLRLoci` had an off-by-one; see [Bugs fixed](#bugs-fixed-vs-the-original).)
- **Resumable runs** — `--checkpoint` saves motif results after the scan; reruns skip the expensive scan and go straight to assembly.
- **Rich reporting** — human-readable summary, TSV statistics (global / per-chromosome / per-motif), and PNG plots out of the box.
- **Graceful interruption** — Ctrl-C outputs whatever batches have completed instead of dropping everything.
- **Original flag compatibility** — drop-in for existing `-i/-x/-y/-o/-g/-b/-m/-a/-f/-c/-t/-n` workflows.

## Quick Start

### Option A — download the prebuilt binary

```bash
# Download from Releases, then:
tar -xzf fastnlr-v1.0.0-linux-x86_64.tar.gz
cd fastnlr-linux-x86_64
./fastnlr --version
./fastnlr -i genome.fasta -p result        # produces result.nlr.txt/.nlr.gff/.nlr.bed/...
```

Verify integrity:

```bash
sha256sum -c fastnlr.sha256
```

### Option B — build from source

```bash
# Requires Rust 1.70+ (developed on 1.94.1)
git clone https://github.com/CropCoder/FastNLR.git
cd NLR-Finder
cargo build --release
# binary: target/release/fastnlr
```

## Usage

```bash
fastnlr -i <input.fasta> [output flags] [options]
```

### Examples

```bash
# 1. Basic loci annotation (uses built-in mot.txt / store.txt)
fastnlr -i genome.fasta -o out.txt -g out.gff -b out.bed

# 2. Prefix-derived subfiles + multithreading
#    -> out.nlr.txt, out.nlr.gff, out.nlr.bed, out.motifs.bed, out.nbarc.fasta
fastnlr -i genome.fasta -p out -t 8

# 3. Full run: report + summary + plots + checkpoint resume
fastnlr -i genome.fasta -p out \
  --stats stats.tsv --summary --plot plots/ --checkpoint ckpt/

# 4. Extract NLR loci sequences (±2000 bp flanking) across ALL contigs
fastnlr -i genome.fasta -p out -f genome.fasta loci.fasta 2000

# 5. Custom motif config override
fastnlr -i genome.fasta -x custom_mot.txt -y custom_store.txt -p out

# 6. Resume from checkpoint after an interrupted run
fastnlr -i genome.fasta -p out --checkpoint ckpt/
```

### Flags

**Input**

| Flag | Description |
|------|-------------|
| `-i <fasta>` | Input genome FASTA (may be gzip-compressed). **Required.** |
| `-x <mot.txt>` | PWM config (optional; default: built-in). |
| `-y <store.txt>` | CDF config (optional; default: built-in). |

**Output** — any combination; omit all to do a scan-only dry run.

| Flag | Description |
|------|-------------|
| `-o <txt>` | NLR loci report (tabular). |
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

## Output Formats

- **`-o` report** — one row per NLR: `seqname, name, domain-class, start, end, strand, motif-list`.
- **`-g` GFF3** — header relabeled for FastNLR with a live system timestamp; `source` column = `FastNLR`.
  ```
  ##gff-version 2
  ##source-version FastNLR V1.0
  ##date 2026-08-20 09:30:01
  ##Type DNA
  ```
- **`-b` BED12** — color-coded blocks (green = complete, orange = partial, red = contains stop codon); reverse-strand blocks reversed.
- **`-a` NB-ARC alignment** — P-loop-anchored multiple alignment with gap-padding for missing motifs.
- **`-f` loci fasta** — extracted ±flanking sequence, reverse-complemented on the reverse strand, 100 bp per line.

## Architecture

FastNLR is a layered Cargo workspace — each crate has a single responsibility and the dependency graph flows bottom-up:

```
nlr-core     domain model + hardcoded rule tables (rank/class/color/seed/signature/consensus)
nlr-config   parse mot.txt (PWM) / store.txt (CDF); built-in embed support
nlr-seq      six-frame translation, reverse complement, codon table, FASTA reader, chopper
nlr-scan     sliding-window scoring + non-overlap arbitration + signature pre-filter (SIMD)
nlr-assemble findSeeds -> mergeSeeds -> elongate three-step assembly
nlr-output   txt/GFF/BED/motifBED/alignment-fasta/loci-fasta/TSV output
nlr-report   run statistics aggregation
nlr-plot     plotters statistics plots (bundled font backend, no system font deps)
nlr-cli      clap CLI entry + pipeline orchestration (rayon + checkpoint + SIGINT); binary: fastnlr
```

**Pipeline:** `FASTA → chop (overlap) → six-frame translate → scan (PWM+CDF) → signature filter → coordinate map → three-step assembly → multi-format output`.

## Bugs fixed vs. the original

1. `-a` P-loop location logic (`while(isPloop)` → `while(!isPloop)`).
2. `-a` stop-codon / unknown-aa replacement (Java `replaceAll` return value was discarded → now actually replaced with `_`).
3. `-f` unified coordinate clamp on the last contig.
4. `-f` multi-contig extraction — the original Java groups NLRs by contig name and extracts from every chromosome; an earlier Rust port only processed the first contig (`seqs.first()`), dropping all other chromosomes. Now correct.
5. `-f` extraction consistency — Java `writeNLRLoci` had an off-by-one (it dropped the first genome character when reading inline, shifting `-f` extraction +1 bp relative to its own GFF/BED coordinates). FastNLR keeps `-f` consistent with the reported coordinates (biologically correct).
6. `##date` header — set to the current system time (Java was non-deterministic).
7. GFF `##source-version` and `source` column relabeled `FastNLR`.
8. p-value rendering in motif BED / export TSV now matches Java `Double.toString` (scientific notation for very small p-values).

## Correctness

Validated against the original algorithm by cross-running the Java jar and diffing every output format. Unit and integration tests cover coordinate mapping (both strands), six-frame translation tail-offset, the codon table, seed/signature tables, scan hits, three-step assembly, and output formats.

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

## Citation

If FastNLR helps your research, please cite the original method:

> Steuernagel, B. et al. *NLR-Annotator: Mining NLRs across the tree of life.* Plant Physiology, 2020. PMID: 32184345.

And reference this Rust port:

> Jiwen Zhao. FastNLR: a high-performance Rust rewrite of NLR-Annotator. https://github.com/CropCoder/FastNLR

## License

GPL-3.0-only — inherited from the original NLR-Annotator. See [LICENSE](LICENSE).

## Author

**Jiwen Zhao** — https://github.com/CropCoder

Contributions and issue reports are welcome.
