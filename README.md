# PLP1 Mutation Mechanism Mapping for Pelizaeus-Merzbacher Disease

A computational pipeline that classifies pathogenic PLP1 variants by molecular mechanism (duplication / loss-of-function / misfolding), validates the misfolding class with structural and transcriptomic evidence, and matches each mechanism to an appropriate therapeutic strategy, rather than assuming a single one-size-fits-all treatment for PMD.

## Summary of Findings

| Phase | What it did | Key result |
|---|---|---|
| 1. Classification | Pulled and classified ClinVar PLP1 variants | 371 pathogenic variants: 173 loss-of-function, 63 duplication, 85 misfolding-candidate |
| 2. Structural analysis | FoldX ΔΔG stability scan on missense variants | 27/83 confirmed destabilizing; worst: Gly246Trp (ΔΔG = 48.3) |
| 3. Transcriptomics | UPR/ER-stress pathway test in a misfolding mutant (GEO GSE111605) | Pathway significantly upregulated, p = 0.0073 |
| 4. Strategy mapping | Combined all evidence into a mechanism-to-therapy map | Full 371-variant table with matched strategy per mechanism |
| 5. Rescue hypothesis | FoldX PositionScan near the worst mutation | Leu31Gly predicted to relieve ~24 kcal/mol of strain |

Full methodology, tools, and per-phase results: see [`WRITEUP.md`](./WRITEUP.md).

## Repository Structure

```text
plp1-pmd-mechanism-mapping/
├── data/
│   ├── clinvar/               # ClinVar variant tables (raw, parsed, classified)
│   └── structures/            # AlphaFold structure, structural scores
├── scripts/
│   ├── phase1_classification/
│   ├── phase2_structural/
│   ├── phase3_transcriptomics/
│   └── phase4_strategy/
├── results/                   # Final output tables (CSV)
├── figures/                   # Generated plots
├── tools/foldx/               # FoldX outputs (binary excluded, see .gitignore)
├── module-b-therapeutics/     # Extended project: therapeutic prediction & intervention design
├── module-c-validation/       # Extended project: twin validation on held-out data
├── PROJECT_SYNTHESIS.md
└── WRITEUP.md                 # Full writeup with methods, tables, figures
```

## Requirements

- Python 3.10 (conda env `pmd-plp1`): biopython, pandas, requests, numpy, matplotlib, scipy
- R (RStudio): DESeq2, GEOquery, limma, biomaRt
- FoldX 5.1 (academic license; download separately from https://foldxsuite.crg.eu/)

## Key Outputs

- `results/plp1_mechanism_strategy_map.csv`: every classified variant with mechanism and matched therapy
- `data/structures/plp1_structural_disruption_scores.csv`: ΔΔG for all 83 missense mutations
- `results/positionscan_gly246trp_rescue.csv`: rescue mutation scan for the worst variant
- `figures/plp1_mutation_map.png`: 3D structure colored by destabilization severity

## Extended Project Modules

- `module-b-therapeutics/`: therapeutic prediction and intervention design
- `module-c-validation/`: twin validation on held-out data (GSE277705)
- Module A (digital twin) is a separate repo: [plp1-pmd-digital-twin](https://github.com/muneeb-dotcom/plp1-pmd-digital-twin)
