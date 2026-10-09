# Test fixtures

`nlr_proteins.fa` holds three full-length NLR proteins used by `tests/scan.rs`:

| Identifier | UniProt | Protein | Expected architecture |
|---|---|---|---|
| `RPM1_ARATH` | Q39214 | Arabidopsis RPM1 | CC(Rx_N)-NB-ARC-LRR (CNL) |
| `RPS4L_ARATH` | Q9SCX6 | Arabidopsis RPS4 | TIR-NB-ARC-LRR (TNL) |
| `ADR1_ARATH` | Q9FW44 | Arabidopsis ADR1 | CC(RPW8)-NB-ARC (RNL) |

Source: UniProtKB (https://rest.uniprot.org/uniprotkb/<acc>.fasta), retrieved 2026-10-09.
UniProt data is distributed under CC BY 4.0; the sequences are reproduced here unmodified
for testing, with the accession preserved in the FASTA header.

The tests additionally generate **shuffled** versions of these sequences at run time and
assert that they are not called as NLRs (negative control).
