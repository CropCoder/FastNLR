# nlr-domain

> HMMER-based protein domain scanning and NLR architecture classification for FastNLR.

`nlr-domain` powers `fastnlr --protein`: it scans protein sequences for NLR characteristic
domains with a pure-Rust HMMER 3.4 port ([`hmmer-pure-rs`](https://docs.rs/hmmer-pure-rs)),
applies Pfam gathering thresholds, and turns the domain evidence into an NLR architecture call.

The crate has no notion of genome coordinates — protein mode works on sequences, not loci.

## Position in the workspace

| | |
|---|---|
| Layer | 3 |
| Depends on | `nlr-core` (domain categories), `hmmer-pure-rs =0.7.4` |
| Used by | `nlr-cli` (`--protein`) |
| Output | `DomainHit` records and `NlrArchitecture` calls |

## Domain evidence: HMMER only for NB-ARC / TIR / CC

Four curated Pfam models are embedded (see `data/README.md` for provenance and licence):

| Accession | Model | Category | GA (sequence / domain) |
|---|---|---|---|
| PF00931 | NB-ARC | NBARC | 23.5 / 23.5 |
| PF01582 | TIR | TIR | 21.3 / 21.3 |
| PF18052 | Rx_N | CC | 27.7 / 27.7 |
| PF05659 | RPW8 | CC | 30.4 / 30.4 |

Two models cover CC because their specificity differs: Rx_N recovers the CC domain of CNL-type
receptors, RPW8 the CC/RPW8 domain of RNL/helper NLRs (ADR1/NRG1/NRC).

**LRR is deliberately not modelled here.** Seven Pfam LRR models combined reach only 47.5%
sensitivity on the RefPlantNLR gold set, while FastNLR's PWM LRR motifs reach 99.5%; the CLI
therefore calls LRR with the PWM scanner and only the three domain families above with HMMER.

Measured on RefPlantNLR (400 proteins, shuffled-sequence negative control):

| Call | Sensitivity | False positives |
|---|---|---|
| NB-ARC | 398/400 (99.5%) | 0/400 |
| TIR | 89/91 (97.8%) | 0/400 |
| CC (N-terminal) | 246/263 (93.5%) | 0/400 |
| LRR (PWM, CLI level) | 391/394 (99.2%) | – |

## Public API

| Item | Description |
|---|---|
| `HmmLibrary::from_embedded()` | the four built-in models |
| `HmmLibrary::add_file(path)` / `add_text(text, label)` | add user models; same accession replaces, new accessions are appended; returns how many models were loaded |
| `HmmLibrary::scanner()` | builds a **thread-local** scanner (must be called on the worker thread) |
| `Scanner::scan_protein(id, seq)` | returns every reported `DomainHit` (passing and failing GA) |
| `DomainModel::ga()` | Pfam GA thresholds from the model's `cutoff` array |
| `DomainHit` | protein id, model accession/name, category, 1-based `start`/`end`, bitscore, p-value, `passed` |
| `classify(hits, lrr_hits)` | `NlrArchitecture { class, nterm, has_nbarc, has_tir, has_lrr, complete, is_nlr }` |
| `category_for_accession(acc)` | map an accession to a `DomainCategory` (unknown → `NA`) |

`hmmer-pure-rs` is pinned exactly (`=0.7.4`) because its scores feed the golden-score tests.

## Classification rules

- Only hits with `passed = true` (sequence **and** domain score above the model's GA) are used.
- `nterm`: the N-terminal-most TIR/CC domain that starts before the first NB-ARC hit;
  otherwise `NB-only` (NB-ARC present) or `None`.
- `class`: categories in positional order with consecutive duplicates collapsed, plus `LRR`
  appended when PWM LRR hits exist — e.g. `CC-NBARC-LRR`, `TIR-NBARC`, `NBARC`.
- `complete`: NB-ARC **and** LRR (protein mode has no P-loop motif anchor, unlike genome mode).
- `is_nlr`: contains NB-ARC, or contains TIR together with LRR.

## Example

```rust
use nlr_domain::{classify, HmmLibrary};

let library = HmmLibrary::from_embedded();
let mut scanner = library.scanner();          // per-thread
let hits = scanner.scan_protein("RPM1", protein_sequence);
let arch = classify(&hits, 3);                 // 3 PWM LRR hits from the caller
assert_eq!(arch.class, "CC-NBARC-LRR");
```

## Testing

```bash
cargo test -p nlr-domain
```

`tests/scan.rs` scans three real NLR proteins (CNL/TNL/RNL, UniProt) with golden bitscores
locked to ±0.5 bit, checks the classification rule table, and asserts that shuffled versions of
the same sequences are not called as NLRs.

## License

GPL-3.0-only for this crate. The bundled Pfam HMMs are CC0 (see `data/README.md`);
`hmmer-pure-rs` is BSD-3-Clause.

Part of [FastNLR](../../README.md) — maintained by Jiwen Zhao ([@CropCoder](https://github.com/CropCoder)).
