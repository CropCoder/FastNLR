# Curated Pfam HMM library for FastNLR protein mode

FastNLR 1.9.1 embeds 23 Pfam profile HMMs. GFF3 mode uses HMM evidence only;
PWM motifs are not calculated.

## Core architecture models

| Accession | Name | Category |
|---|---|---|
| PF00931 | NB-ARC | NBARC |
| PF01582 | TIR | TIR |
| PF18052 | Rx_N | CC |
| PF05659 | RPW8 | CC |
| PF12061 | R1-type / NB-LRR | CC |

The three CC models are retained separately because Rx_N, RPW8 and R1-type
domains distinguish CNL, RNL/helper and NRC-clade architectures.

## LRR models

PF00560, PF07723, PF07725, PF12799, PF13516, PF13855, PF14580, PF18805,
PF18831, PF18837, PF23598 and PF24969 cover generic LRR repeats.
PF23247 and PF23286 are plant NLR-specific LRR profiles. PF01462, PF08263,
PF01463 and PF16920 cover LRR N- or C-terminal regions.

The LRR models are scanned only after a translated protein has passed GA for
NB-ARC or TIR. A passed LRR is C-terminal evidence only when its protein start
is after the end of the first NB-ARC block. LRR hits before NB-ARC remain in
the domain table but are excluded from the architecture and cannot make a
protein complete.

## Exact model snapshot

Source: InterPro/Pfam HMM API, retrieved 2026-11-11:
`https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<ACC>/?annotation=hmm`.
The five core models retain their validated 1.9.0 versions to avoid a
regression observed when refreshing to PF00931.30 (NB-ARC coverage decreased
from 398/400 to 395/400). `fetch_hmms.py` keeps these pinned core files unless
`FASTNLR_REFRESH_CORE=1` is explicitly set.

| Accession | Version | Name | GA sequence/domain | SHA256 |
|---|---|---|---|---|
| PF00560 | PF00560.41 | LRR_1 | 20.6 / 9.3 | `46e3385b79ab2bc7ab3356e8faf8de29f4b7310e5a601411e8598743d902a65b` |
| PF00931 | PF00931.29 | NB-ARC | 23.5 / 23.5 | `c83b56aa0ccbc701372eb8ca434df9d527287f2824ea874fb8370c1d465a2326` |
| PF01462 | PF01462.24 | LRRNT | 25.8 / 24.4 | `66af4d215449714eb18e92bba55a58582919452e3f9373e4414409bfc328ad55` |
| PF01463 | PF01463.31 | LRRCT | 20.0 / 6.6 | `62a74be7679c644ae9ee6b1a35bcf37acb1942a10be9b044b2de458de7adac13` |
| PF01582 | PF01582.26 | TIR | 21.3 / 21.3 | `3253c338c8c5d4a0b4b81682a35f2dc73731f2de35ca79448419e1a7272f8f1b` |
| PF05659 | PF05659.17 | RPW8 | 30.4 / 30.4 | `7036821a4df881e42d5bad641abf44060d3bff626e9493e71cbf27163da97561` |
| PF07723 | PF07723.20 | LRR_2 | 20.9 / 9.5 | `f7b868d4fde686728b0152e9d5b4be39719c768d9d351c152fbddcc3a3d50aff` |
| PF07725 | PF07725.19 | LRR_3 | 20.2 / 20.2 | `f2e7b979da03a7d817266306728c2e22b9fcb669a06ce1a1141b0a1e2d52b390` |
| PF08263 | PF08263.19 | LRRNT_2 | 20.7 / 3.1 | `00d4bbd467d20f5af989dda63acaf8e556e6b86aa0420732cf788461f539853b` |
| PF12061 | PF12061.14 | NB-LRR | 28.7 / 28.7 | `0d3853d3f930b72dc93dbdd3f40ab65165887302a80b7b959a9454093a3fe361` |
| PF12799 | PF12799.15 | LRR_4 | 27.0 / 27.0 | `7ce5af95250c5acc9caf8becb09d64a6bf4f544b99c82cbfe8aa811f8cc4477d` |
| PF13516 | PF13516.14 | LRR_6 | 23.4 / 8.7 | `02579e7eda65360c4a5ecebe6857e7760a306fbb65103c8dc1d04221fe468111` |
| PF13855 | PF13855.14 | LRR_8 | 27.0 / 27.0 | `8eaa390fa4246f1162b673c2477d55bf0802823dad390b66f34ef8ccbd9438ff` |
| PF14580 | PF14580.13 | LRR_9 | 30.5 / 30.5 | `c1713c47de6c11354c6eddc8a507af260ac79372677263f657c99b7575ef9295` |
| PF16920 | PF16920.12 | LRRCT_2 | 23.8 / 23.8 | `ee261cb039b6290d10932fa6a3f9f0ba864b4ea8b0bec653735f2a6f389f35c1` |
| PF18052 | PF18052.7 | Rx_N | 27.7 / 27.7 | `7c7d423fd1c679a28c55b527eaa2fcf3b777cf8af5a750b413fd3802e3104b22` |
| PF18805 | PF18805.9 | LRR_10 | 30.0 / 20.0 | `a62b08f029c08f9b42e7a10031a5705869c8aafee26fb7cec3687aee0f491722` |
| PF18831 | PF18831.9 | LRR_11 | 29.0 / 20.0 | `8a540786b4ee2b311112155e4267971bd829a99db25483094474103019ebf267` |
| PF18837 | PF18837.9 | LRR_12 | 34.0 / 20.0 | `288a66f30b40a849dab46b6f0ec219d4c716590f84d98a219968816480caa8f3` |
| PF23247 | PF23247.3 | LRR_RPS2 | 27.9 / 27.9 | `76396d1231a362ce37892fe3056255f912e06acc27cb9e804d10dafe93f0a98c` |
| PF23286 | PF23286.1 | LRR_13 | 30.1 / 30.1 | `d00c1689b945154e2656415637b7db52a9cd2cc1f645b9e1f603c056c1d11bf4` |
| PF23598 | PF23598.2 | LRR_14 | 40.0 / 40.0 | `bfe9bc77d63b8ba769ddfc8bf18bdf3feb02d0e10782690b4f8d90b59facd082` |
| PF24969 | PF24969.2 | LRR_15 | 24.0 / 24.0 | `c88929409a02a12622193d6ec10a78ab121a601dbe8f6097fa7b077fa7091d4c` |

## Validation and provenance

The library was checked against the 400 RefPlantNLR proteins. The five core
models detected NB-ARC in 398/400 proteins; the 18 LRR models jointly detected
at least one LRR in 349/394 proteins annotated with an LRR in the RefPlantNLR
structure table. No model passed GA on the 400 deterministically shuffled
negative controls. The ordered classifier marked 364/394 proteins with
NB-ARC plus LRR in the truth table as `complete` (92.4%).

The formal validation is implemented in
`benchmark/scripts/S42-Validate-NLR-Domain-HMM-Library-V1.R`.

Methodological references:

- Kourelis J, Sakai T, Adachi H, Kamoun S. RefPlantNLR is a comprehensive
  collection of experimentally validated plant disease resistance proteins
  from the NLR family. PLoS Biology (2021). PMID 34669691.
- Barragan AC, Weigel D. Plant NLR diversity: the known unknowns of
  pan-NLRomes. Plant Cell (2021). PMID 33793812.
- Kobe B, Deisenhofer J. The leucine-rich repeat: a versatile binding motif.
  Trends in Biochemical Sciences (1994). PMID 7817399.
- Steuernagel B, et al. The NLR-Annotator Tool Enables Annotation of the
  Intracellular Immune Receptor Repertoire. Plant Physiology (2020).
  PMID 32184345.
- Paysan-Lafosse T, et al. The Pfam protein families database: embracing
  AI/ML. Nucleic Acids Research (2025). PMID 39540428.
- Blum M, et al. InterPro: the protein sequence classification resource in
  2025. Nucleic Acids Research (2025). PMID 39565202.
- Castel B, et al. Diverse NLR immune receptors activate defence via the
  RPW8-NLR NRG1. New Phytologist (2019). PMID 30582759.
- Contreras MP, et al. NLR receptors in plant immunity: making sense of the
  alphabet soup. EMBO Reports (2023). PMID 37602936.

Pfam/InterPro data is distributed under CC0 1.0. The models are embedded
unmodified; only the gzip wrapper is removed.

## Regenerating

```bash
python3 fetch_hmms.py
```

After regenerating, run:

```bash
cargo test -p nlr-domain --locked
Rscript benchmark/scripts/S42-Validate-NLR-Domain-HMM-Library-V1.R --rerun
```
