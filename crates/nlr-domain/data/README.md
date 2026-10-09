# Curated domain models for FastNLR protein mode

These four Pfam profile HMMs are embedded into the `nlr-domain` crate and used by
`fastnlr --protein` to locate NLR characteristic domains.

| File | Pfam accession | Name | Category | GA (sequence / domain) |
|---|---|---|---|---|
| `PF00931.hmm` | PF00931.29 | NB-ARC | NBARC | 23.5 / 23.5 |
| `PF01582.hmm` | PF01582.26 | TIR | TIR | 21.3 / 21.3 |
| `PF18052.hmm` | PF18052.7 | Rx_N | CC | 27.7 / 27.7 |
| `PF05659.hmm` | PF05659.17 | RPW8 | CC | 30.4 / 30.4 |

## Why these models

- **PF00931 (NB-ARC)** and **PF01582 (TIR)** are the diagnostic domains of the two major
  NLR branches; both are used with their Pfam gathering thresholds.
- **PF18052 (Rx_N)** covers the CC domain of CNL-type receptors, **PF05659 (RPW8)** covers the
  CC/RPW8 domain of RNL/helper NLRs (ADR1/NRG1/NRC). They are kept as separate models because
  their specificity differs: on the RefPlantNLR gold set Rx_N recovers 249/263 CC-type
  proteins while RPW8 only fires on the RNL subset.
- **LRR is deliberately not modelled here.** Seven Pfam LRR models combined reach only 47.5%
  sensitivity on RefPlantNLR NLRs, whereas FastNLR's existing PWM motifs reach 99.5%. Protein
  mode therefore calls LRR with the PWM motifs (motif_9/motif_11) and only the domains above
  with HMMER.

## Provenance and licence

- Source: InterPro REST API, `https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<acc>/?annotation=hmm`
  (retrieved 2026-10-09; the API returns a gzip-compressed HMMER3 ASCII file).
- Pfam/InterPro data is distributed under **CC0 1.0** (public domain dedication); the HMM
  files may be redistributed. Attribution is given here regardless.
- The models are used unmodified; only the gzip wrapper is removed.

## Regenerating

```bash
python3 fetch_hmms.py            # writes into this directory
```

After regenerating, re-run `cargo test -p nlr-domain`: the golden-score assertions
(±0.5 bit) will fail if a model changed, which is the intended signal to re-validate.
