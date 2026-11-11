# nlr-domain

> HMMER-based protein domain scanning and NLR architecture classification for FastNLR.

`nlr-domain` powers FastNLR's GFF3-guided protein workflow: it scans translated gene proteins for NLR characteristic
domains with a pure-Rust HMMER 3.4 port ([`hmmer-pure-rs`](https://docs.rs/hmmer-pure-rs)),
applies Pfam gathering thresholds, and turns the domain evidence into an NLR architecture call.

The crate works on protein sequences; `nlr-cli` maps protein coordinates through the selected
CDS model back to genomic coordinates.

## Position in the workspace

| | |
|---|---|
| Layer | 3 |
| Depends on | `nlr-core` (domain categories), `hmmer-pure-rs =0.7.4` |
| Used by | `nlr-cli` (`--gff`) |
| Output | `DomainHit` records and `NlrArchitecture` calls |

## Domain evidence: HMMER only

Twenty-three Pfam models are embedded (see `data/README.md` for provenance,
versions, GA thresholds and SHA256):

| Accession | Model | Category | GA (sequence / domain) |
|---|---|---|---|
| PF00931 | NB-ARC | NBARC | 23.5 / 23.5 |
| PF01582 | TIR | TIR | 21.3 / 21.3 |
| PF18052 | Rx_N | CC | 27.7 / 27.7 |
| PF05659 | RPW8 | CC | 30.4 / 30.4 |
| PF12061 | R1-type (Pfam short name `NB-LRR`) | CC | 28.7 / 28.7 |
| PF00560, PF07723, PF07725, PF12799, PF13516, PF13855, PF14580, PF18805, PF18831, PF18837, PF23598, PF24969 | generic LRR repeats | LRR | model-specific GA |
| PF23247, PF23286 | plant NLR-specific LRR profiles | LRR | model-specific GA |
| PF01462, PF08263, PF01463, PF16920 | LRR N/C-terminal regions | LRR | model-specific GA |

Three models cover the CC family because their specificity differs: Rx_N recovers the CC domain
of most CNL-type receptors, RPW8 the CC/RPW8 domain of RNL/helper NLRs (ADR1/NRG1), and PF12061
the R1-type N-terminal domain of the NRC-clade CNLs (RefPlantNLR `Structure` letter `B`;
NbPRFa/NbPRFb/Prf/R1).

LRR models are scanned only after NB-ARC or TIR passes GA. This preserves
the common non-NLR proteome fast path while giving NLR candidates explicit
position-aware LRR evidence.

Measured on RefPlantNLR (400 proteins, shuffled-sequence negative control):

| Call | Sensitivity | False positives |
|---|---|---|
| NB-ARC | 398/400 (99.5%) | 0/400 |
| TIR | 89/91 (97.8%) | 0/400 |
| CC (N-terminal) | 246/263 (93.5%) | 0/400 |
| LRR (18 HMM union) | 349/394 (88.6%) | 0/400 shuffled |

## Public API

| Item | Description |
|---|---|
| `HmmLibrary::from_embedded()` | all 23 built-in models |
| `HmmLibrary::add_file(path)` / `add_text(text, label)` | add user models; same accession replaces, new accessions are appended; returns how many models were loaded |
| `HmmLibrary::core_scanner()` | thread-local scanner for five core models |
| `HmmLibrary::protein_scanner()` | thread-local two-stage core/LRR scanner |
| `Scanner::scan_protein(id, seq)` | returns every reported `DomainHit` (passing and failing GA) |
| `DomainModel::ga()` | Pfam GA thresholds from the model's `cutoff` array |
| `DomainHit` | protein id, model accession/name, category, 1-based `start`/`end`, bitscore, p-value, `passed` |
| `classify(hits)` | `NlrArchitecture { class, family, nterm, has_nbarc, has_tir, has_lrr_any, has_cterm_lrr, complete, is_nlr }` |
| `category_for_accession(acc)` | map an accession to a `DomainCategory` (unknown → `NA`) |

`hmmer-pure-rs` is pinned exactly (`=0.7.4`) because its scores feed the golden-score tests.

## Classification rules

- Only hits with `passed = true` (sequence **and** domain score above the model's GA) are used.
- `nterm`: the N-terminal-most TIR/CC domain that starts before the first NB-ARC hit;
  otherwise `NB-only` (NB-ARC present) or `None`.
- `class`: categories in positional order with consecutive duplicates collapsed.
  LRR before NB-ARC is ignored; LRR after NB-ARC is retained at its real position.
- `has_cterm_lrr`: at least one passed LRR starts after the first NB-ARC block.
- `complete`: NB-ARC and C-terminal LRR are both present.
- `family`: biological lineage inferred from the N-terminal domain and NB-ARC.
- `is_nlr`: NB-ARC is present, or TIR is present together with C-terminal LRR.

## Example

```rust
use nlr_domain::{classify, HmmLibrary};

let library = HmmLibrary::from_embedded();
let mut scanner = library.protein_scanner(); // per-thread
let hits = scanner.scan_protein("RPM1", protein_sequence);
let arch = classify(&hits);
assert_eq!(arch.class, "CC-NBARC-LRR");
assert!(arch.complete);
```

## Testing

```bash
cargo test -p nlr-domain
```

`tests/scan.rs` scans four real NLR proteins (CNL/TNL/RNL plus an NRC-clade CNL, all from
UniProt) with golden bitscores locked to ±0.5 bit, checks the classification rule table, and
asserts that shuffled versions of the same sequences are not called as NLRs. The
`NlrArchitecture.family` field distinguishes RNL from CNL even when both share a
`CC-NBARC-LRR` architecture: RPW8(PF05659) gives `RNL`, while Rx_N/R1-type CC gives `CNL`.

## License

GPL-3.0-only for this crate. The bundled Pfam HMMs are CC0 (see `data/README.md`);
`hmmer-pure-rs` is BSD-3-Clause.

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
