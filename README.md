\# PLP1 Mutation Mechanism Mapping for Pelizaeus-Merzbacher Disease



A computational pipeline classifying pathogenic PLP1 variants by molecular

mechanism (duplication / loss-of-function / misfolding), validating the

misfolding class with structural and transcriptomic evidence, and matching

each mechanism to an appropriate therapeutic strategy — rather than

assuming a single one-size-fits-all treatment for PMD.



\## Summary of Findings



| Phase | What it did | Key result |

|---|---|---|

| 1. Classification | Pulled \& classified ClinVar PLP1 variants | 371 pathogenic variants → 173 loss-of-function, 63 duplication, 85 misfolding-candidate |

| 2. Structural analysis | FoldX ΔΔG stability scan on missense variants | 27/83 confirmed destabilizing; worst: Gly246Trp (ΔΔG=48.3) |

| 3. Transcriptomics | UPR/ER-stress pathway test in a misfolding mutant (GEO GSE111605) | Pathway significantly upregulated, p=0.0073 |

| 4. Strategy mapping | Combined all evidence into mechanism→therapy map | Full 371-variant table with matched strategy per mechanism |

| 5. Rescue hypothesis | FoldX PositionScan near worst mutation | Leu31Gly predicted to relieve \~24 kcal/mol of strain |



Full methodology, tools, and per-phase results: see \[`WRITEUP.md`](./WRITEUP.md)



\## Repository Structure

plp1-pmd-project/

├── data/

│ ├── clinvar/ # ClinVar variant tables (raw, parsed, classified)

│ └── structures/ # AlphaFold structure, structural scores

├── scripts/

│ ├── phase1\_classification/

│ ├── phase2\_structural/

│ ├── phase3\_transcriptomics/

│ └── phase4\_strategy/

├── results/ # Final output tables (CSV)

├── figures/ # Generated plots

├── tools/foldx/ # FoldX outputs (binary excluded, see .gitignore)

└── WRITEUP.md # Full writeup with methods, tables, figures



\## Requirements



\- Python 3.10 (conda env `pmd-plp1`): biopython, pandas, requests, numpy, matplotlib, scipy

\- R (RStudio): DESeq2, GEOquery, limma, biomaRt

\- FoldX 5.1 (academic license, must be downloaded separately from https://foldxsuite.crg.eu/)



\## Key Outputs



\- `results/plp1\_mechanism\_strategy\_map.csv` — every classified variant with mechanism + matched therapy

\- `data/structures/plp1\_structural\_disruption\_scores.csv` — ΔΔG for all 83 missense mutations

\- `results/positionscan\_gly246trp\_rescue.csv` — rescue mutation scan for the worst variant

\- `figures/plp1\_mutation\_map.png` — 3D structure colored by destabilization severity

