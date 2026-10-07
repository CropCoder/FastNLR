# nlr-core

> Domain model and static biology rule tables for the FastNLR workspace.

`nlr-core` is the foundation of FastNLR. It has **no external dependencies and performs no I/O**:
it defines how a motif hit and an NLR locus are represented in memory, how protein coordinates
are mapped back onto genome coordinates, and which hard-coded rules (rank, domain category,
colour, seed combinations, signatures, consensus sequences) drive seeding, assembly and
reporting.

Every other crate reads its rules from here, so changing a table in this crate propagates
consistently through scanning, assembly, every output format and every plot.

## Position in the workspace

| | |
|---|---|
| Layer | 1 (bottom) |
| Depends on | nothing |
| Used by | `nlr-config`, `nlr-scan`, `nlr-assemble`, `nlr-output`, `nlr-report`, `nlr-cli` |
| I/O | none — pure data structures and rule evaluation, fully unit-testable |

## Modules

| Module | Responsibility |
|---|---|
| `strand` | `Strand` enum (`Forward` / `Reverse`) with `symbol()` for output formats |
| `motif` | `Motif` — a single motif hit, carrying both protein and genome coordinates |
| `motif_list` | `MotifList` — a named, ordered motif collection representing one locus |
| `signature_def` | Static rule tables and `AnnotatorSignatureDefinition` |

## Core concepts

### Coordinate model

A `Motif` starts life on the protein side (`new_protein`) and is promoted to genome coordinates
by `set_dna`, which is called once the fragment offset, reading frame and strand are known:

```
forward:  dna_start = (position - 1) * 3 + frame + offset
reverse:  dna_start = offset + fragment_length - ((position + len - 1) * 3 + frame)
```

`position` is 1-based and `len` is the motif length (the matched protein string). On **both**
strands `dna_start` is the leftmost coordinate and `dna_end` the rightmost, so `dna_start < dna_end`
always holds and downstream code never has to special-case the strand.

### Asymmetric ordering

`Motif`'s `Ord` implementation is deliberately asymmetric, because a locus is walked in
translation direction rather than in genome order:

- same contig, forward strand — ascending `dna_start`;
- same contig, reverse strand — descending `dna_start`;
- different contigs — ordered by `protein_sequence_id`;
- DNA parameters not yet set — ordered by `protein_sequence_id`, then `position`.

`MotifList` keeps its motifs sorted at all times, so `first_motif()` and `last_motif()` always
mean "first and last in translation direction".

### Locus semantics

| Predicate | Meaning |
|---|---|
| `is_complete_nlr` | contains the P-loop motif (`motif_1`) **and** at least one LRR motif |
| `has_stop_codon` | any constituent motif contains `*` in its translated sequence |
| `can_be_merged_with` | same strand, one of the four `dna_start` endpoint combinations within the threshold, and rank stays monotonically non-decreasing after merge (LRR motifs may tie) |
| `span` | min/max over the first and last motif's `dna_start`/`dna_end` |

### Rule tables

All tables are compile-time constants, indexed by `MotifId` (`u8`). Index `0` is an unused
placeholder, ids `1..=20` are built-in, and ids `21..=28` are reserved for the optional
RNL/TIR extension library shipped in `data/motif_library_rnl_tir/`.

| Constant | Size | Purpose |
|---|---|---|
| `PLOOP_MOTIF` | 1 | id of the P-loop motif (`motif_1`), the NB-ARC anchor used by the alignment output |
| `BUILTIN_MOTIF_COUNT` | 20 | ids above this value are library-specific and default to `DomainCategory::Na` |
| `RANKS` | 29 | translation-direction ordering used for seeding and elongation |
| `CATEGORIES` | 29 | `NBARC` / `LRR` / `TIR` / `CC` / `LINKER` / `NA` |
| `RGB_COLORS` | 29 | per-motif colour emitted into BED `itemRgb` |
| `SEED_COMBINATIONS` | 13 | motif-id sequences that trigger seed creation |
| `SIGNATURES` | 18 | local cluster pre-filter patterns |
| `CONSENSUS_SEQUENCES` | 29 | gap lengths used by the NB-ARC multiple alignment |

## Public API

### `Motif`

| Item | Description |
|---|---|
| `new_protein(id, protein_sequence_id, position, protein_sequence, pvalue)` | constructs a protein-side hit |
| `set_dna(dna_sequence_id, offset, fragment_length, frame, strand)` | maps to genome coordinates and sets `dna_parameters_set` |
| `has_stop()` | whether the matched protein string contains `*` |
| `export_string()` | serialises 11 TSV fields (id, protein id, position, sequence, p-value, score, contig, start, end, strand, frame) |
| `from_export_line(line)` | inverse of `export_string`, used by checkpoint import |
| `format_double_java(d)` | formats `f64` in plain notation for `1e-3 <= abs(d) < 1e7`, otherwise upper-case scientific notation |

### `MotifList`

| Item | Description |
|---|---|
| `new(name, motifs)` / `sort()` | construction and re-sorting |
| `add_motif` / `add_motifs` | append and re-sort |
| `first_motif` / `last_motif` / `strand` / `is_forward` / `sequence_name` | accessors in translation direction |
| `is_complete_nlr(def)` / `has_stop_codon()` | locus classification |
| `can_be_merged_with(other, distance, def)` | merge predicate used by `nlr-assemble` |
| `remove_redundant_motifs()` | deduplicate by `(id, dna_start, dna_end)` |
| `motif_list_string()` / `domain_string(def)` / `span()` | reporting helpers |

### `AnnotatorSignatureDefinition`

Built from the static tables and extended through a consuming builder:

```rust
let def = AnnotatorSignatureDefinition::new()
    .with_extra_categories(extra_categories)
    .with_extra_rules(extra_seeds, extra_signatures)
    .with_flexible_seed(FlexibleSeedConfig { min_nbarc: 3, require_ploop: false });
```

Queries: `rank`, `category`, `is_lrr`, `is_nbarc`, `is_ploop`, `color_rgb`, `consensus`,
`nbarc_motif_order`, `is_seed`, `is_seed_relaxed`, `is_seed_flexible`, `has_signature`,
`has_signature_flexible`, `domain_string`.

Three seeding strategies are supported, in increasing order of recall for divergent loci:

1. `is_seed` — matches `SEED_COMBINATIONS` (plus CLI-declared extras). Built-in combinations that
   only reference ids `<= BUILTIN_MOTIF_COUNT` require an **exact** match of the accumulated id
   sequence; combinations that reference an external-library id (`> 20`) match as a **contiguous
   subsequence** instead, because extended-library hits are additionally prefixed by their own
   library-specific motifs;
2. `is_seed_relaxed(seq, min_nbarc)` — a consecutive run containing at least `min_nbarc` NB-ARC motifs;
3. `is_seed_flexible(seq, cfg)` — category/rank-based: `min_nbarc` NB-ARC motifs appearing in rank order, optionally requiring the P-loop.

Two further queries drive filtering and reporting:

- `has_signature(ids)` — whether any signature occurs as a **contiguous subsequence**;
- `domain_string(ids)` — collapses consecutive equal categories, skips `NA` and `LINKER`, and
  joins the remainder with `-` (for example `CC-NBARC-LRR`). This is the string written to the
  `nlrClass` attribute of the GFF and to the third column of the text report.

## Example

```rust
use nlr_core::{AnnotatorSignatureDefinition, Motif, Strand};

let def = AnnotatorSignatureDefinition::new();

let mut hit = Motif::new_protein(1, "chr1_0_frame+0".to_string(), 42, "GxxxxGKS".to_string(), 1e-8);
hit.set_dna("chr1".to_string(), 0, 20_000, 0, Strand::Forward);

assert!(def.is_ploop(hit.id));
assert!(hit.dna_start < hit.dna_end);
assert_eq!(def.motif_id_str(21), "motif_21");
```

## Invariants

- `dna_start < dna_end` for every motif with `dna_parameters_set`.
- `MotifList::motifs` is always sorted in translation direction after construction or mutation.
- Rule tables are never mutated at runtime; external libraries extend behaviour through
  `extra_categories` / `extra_seeds` / `extra_signatures` instead of patching constants.
- Motif ids are `u8`; ids above `BUILTIN_MOTIF_COUNT` are legal but require an explicit category.

## Testing

```bash
cargo test -p nlr-core
```

`tests/core.rs` covers coordinate mapping on both strands, rule-table lookups, seed/signature
predicates and the flexible seeding path; `tests/regression.rs` pins behaviour that must not
change between releases.

## License

GPL-3.0-only. See the workspace [LICENSE](../../LICENSE).

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
