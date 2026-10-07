# nlr-seq

> FASTA ingestion, sequence chopping, six-frame translation and reverse complement.

`nlr-seq` is the sequence layer of FastNLR. It streams a genome FASTA in overlapping fragments,
translates every fragment in all six reading frames, and exposes the small helpers (codon
translation, reverse complement, DNA detection) that the rest of the pipeline needs.

It is a standalone crate: it depends only on `flate2` and `memmap2`, and knows nothing about
motifs or NLRs.

## Position in the workspace

| | |
|---|---|
| Layer | 2 |
| Depends on | `flate2`, `memmap2` |
| Used by | `nlr-scan`, `nlr-cli` |
| I/O | streaming reads; memory-mapped fast path for uncompressed input |

## Modules

| Module | Responsibility |
|---|---|
| `codon` | The 64-entry standard genetic code and `translate_triplet` |
| `translate` | `BioSequence` with six-frame translation and reverse complement |
| `fasta` | `FastaReader` (streaming, gzip-aware) and `read_all` |
| `chopper` | `SequenceChopper` — overlapping fragment stream |

## Public API

### `codon`

| Item | Description |
|---|---|
| `translate_triplet(&[u8]) -> char` | standard genetic code; `TAA`/`TAG`/`TGA` → `*`; any triplet containing a non-`ACGT` base → `X` |

### `translate::BioSequence`

| Item | Description |
|---|---|
| `new(identifier, sequence)` | constructs an entry; `identifier` is the first whitespace-delimited header token |
| `len()` / `is_empty()` | sequence length in bytes |
| `reverse_complement()` | case-preserving complement, non-standard bases passed through |
| `translate2protein()` | six frames as `{id}_frame+0..+2` and `{id}_frame-0..-2` |
| `is_dna()` | heuristic: `A`/`T`/`G`/`C`/`N` make up more than 50% of the first 500 characters |
| `fasta_string()` | re-emits the entry as FASTA (with description) |

Six-frame translation reproduces the reference implementation's frame-specific tail semantics
through a small `end_offset` table derived from `len % 3`, so the last partial codon is handled
identically on both strands.

### `fasta`

| Item | Description |
|---|---|
| `FastaReader::from_file(path)` | detects gzip by magic bytes (`1f 8b`) and wraps the file in a multi-member gzip decoder when needed |
| `FastaReader::new(reader)` | reader-agnostic constructor for tests and in-memory input |
| `read_entry()` | yields the next `BioSequence`, `None` at EOF; the description keeps the remainder of the header |
| `read_all(path)` | loads every entry into memory — convenient, but peak memory scales with the genome |

### `chopper::SequenceChopper`

| Item | Description |
|---|---|
| `from_file(path, fragment_length, overlap)` | sniffs gzip; plain FASTA takes the zero-copy `mmap` path, gzip takes the streaming path |
| `from_mmap_file` / `new` | explicit constructors for the two paths |
| `next_sequence()` | yields fragments named `{contig}_{offset}`, uppercased, `None` at EOF |

Fragments advance by `fragment_length - overlap`, so consecutive fragments share `overlap` bases.
The overlap exists to guarantee that every motif is fully contained in at least one fragment,
including motifs that straddle a fragment boundary; it is deliberately independent of how loci
are assembled later.

## Example

```rust
use nlr_seq::{BioSequence, SequenceChopper};

let seq = BioSequence::new("chr1", "ATGAAATAA");
let frames = seq.translate2protein();
assert_eq!(frames.len(), 6);
assert_eq!(frames[0].identifier, "chr1_frame+0");
assert!(frames[0].sequence.starts_with("MK"));
```

## Design notes

- **Two ingestion paths, one behaviour.** The `mmap` path avoids a copy per fragment, while the
  gzip path decompresses on the fly; both produce identical fragment identifiers and sequences.
- **Fragment identifiers are the coordinate contract.** `{contig}_{offset}` is parsed by
  `nlr-cli` to recover the contig name and the fragment offset when mapping motif hits back to
  genome coordinates.
- **`read_all` is a convenience, not the streaming path.** Callers that process whole genomes
  should prefer `SequenceChopper`; `read_all` exists for small inputs and for bounded uses.

## Testing

```bash
cargo test -p nlr-seq
```

Unit tests in `src/codon.rs`, `src/translate.rs` and `src/chopper.rs` cover the codon table,
reverse complement, six-frame offsets and fragment stepping; `tests/mmap.rs` and
`tests/mmap_edge.rs` cover the memory-mapped path, including malformed headers and EOF handling.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
