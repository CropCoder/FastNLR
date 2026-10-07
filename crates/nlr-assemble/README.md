# nlr-assemble

> Three-step assembly of motif hits into NLR loci.

`nlr-assemble` turns a flat, coordinate-sorted list of motif hits on one contig into named NLR
loci. It implements the `findSeeds → mergeSeeds → elongate` algorithm: find the motif
combinations that are diagnostic for an NLR, merge compatible seeds that describe the same
gene, then extend both ends along the translation direction while the domain order stays
consistent.

The crate depends only on `nlr-core`; it performs no I/O and is fully deterministic.

## Position in the workspace

| | |
|---|---|
| Layer | 3 |
| Depends on | `nlr-core` |
| Used by | `nlr-cli` |
| Input | all motif hits of a single contig |
| Output | named `MotifList` values, one per locus |

## Algorithm

### Step 1 — `find_seeds`

Walk the sorted hits and accumulate a candidate motif-id sequence while consecutive hits stay
within `distance_within_motif_combination`. After each extension the sequence is tested against
the seeding strategy:

| Mode | Trigger |
|---|---|
| exact (default) | the id sequence matches one of `SEED_COMBINATIONS` |
| relaxed (`--relaxed-seed n`) | a consecutive run containing at least `n` NB-ARC motifs |
| flexible (`--flexible-seed n`) | at least `n` NB-ARC motifs appear in rank order (optionally requiring the P-loop) |

The first match closes the seed and the walk resumes from the next hit.

### Step 2 — `merge_seeds`

Seeds are sorted and then merged pairwise while `MotifList::can_be_merged_with` holds: same
strand, at least one of the four `dna_start` endpoint combinations within
`distance_between_motif_combinations`, and a monotonically non-decreasing rank order once the
combined list is deduplicated and sorted (LRR motifs may tie). Merging is what collapses
multiple seeds belonging to the same gene into a single locus.

### Step 3 — `elongate`

Each merged seed grows outward from its first and last motif while the next hit is within
`distance_for_elongating` and continues the domain order:

- forward direction: `rank(candidate) < rank(current first motif)`;
- reverse direction: `rank(candidate) > rank(current last motif)`, or equal rank for an LRR motif.

A motif containing a stop codon is only admitted when it is an NB-ARC motif. Every motif
consumed by a locus is recorded in a `used` set keyed by `(id, dna_start)`, and a seed whose
motifs have already been used is skipped, so a hit belongs to at most one locus.

### Naming

Candidate loci are sorted in the crate-wide translation-direction order before naming, and each
receives `{seq_id}_nlr{N}` with `N` starting at 1. Naming is therefore stable across runs and
independent of thread scheduling.

## Public API

| Item | Description |
|---|---|
| `AssembleParams` | tuning parameters, see below |
| `assemble(seq_id, motifs, params, def) -> Vec<MotifList>` | runs all three steps and returns the named loci |

### `AssembleParams`

| Field | Default | Meaning |
|---|---|---|
| `distance_within_motif_combination` | 500 | maximum gap between adjacent motifs inside one seed |
| `distance_for_elongating` | 2 500 | maximum gap while extending a seed |
| `distance_between_motif_combinations` | 10 000 | maximum gap for merging two seeds |
| `relaxed_seed_min_nbarc` | `None` | enables relaxed seeding with the given NB-ARC count |

## Example

```rust
use nlr_assemble::{assemble, AssembleParams};
use nlr_core::AnnotatorSignatureDefinition;

let def = AnnotatorSignatureDefinition::new();
let params = AssembleParams::default();

let loci = assemble("chr1", motifs, &params, &def);
for locus in &loci {
    println!("{} {}..{}", locus.name, locus.span().0, locus.span().1);
}
```

## Guarantees

- **Determinism.** The input is sorted on entry and every intermediate vector is sorted before
  iteration, so identical input always yields identical locus names and contents.
- **No double assignment.** The `(id, dna_start)` used-set prevents a motif from being claimed by
  two loci.
- **No implicit library assumptions.** All rank/category/seed rules come from
  `AnnotatorSignatureDefinition`, so external motif libraries change behaviour through the
  definition rather than through this crate.

## Performance notes

Assembly is linear-ish in the number of hits per contig, with a quadratic component in
`merge_seeds` (seeds are merged from the front of a sorted vector). Contigs are independent, and
`nlr-cli` assembles them in parallel.

## Testing

```bash
cargo test -p nlr-assemble
```

`tests/assemble.rs` covers seeding on exact combinations, seed merging, elongation in both
directions, stop-codon handling and locus naming.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
