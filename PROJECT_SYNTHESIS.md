# PLP1/PMD Extended Project — Full Synthesis

## 1. PLP1 Mechanistic Heterogeneity (Part 1)
371 pathogenic PLP1 variants classified into 3 mechanism classes: 173 loss-of-function, 63 duplication, 85 misfolding candidates. FoldX ddG scan confirmed 27 as structurally destabilizing (top: Gly246Trp, ddG=48.3). UPR pathway significantly upregulated in a misfolding-mutant model (Wilcoxon p=0.0073). FoldX PositionScan identified Leu31Gly as a candidate rescue mutation (-24.3 kcal/mol).

## 2. Digital Twin (Module A)
Built a mouse oligodendrocyte-lineage reference atlas (Marques et al. 2016, 4,927 QC'd cells), validated against canonical markers, with a 6-module cell-state framework (myelin_output, er_stress, opc_identity, apoptosis, lipid_synth) and a maturation pseudotime axis. Modules validated on independent human data (Jakel et al. 2019). Cross-mechanism expression data confirmed unavailable for LOF/duplication classes -- twin honestly reports 'no_measured_data' for those rather than extrapolating. Full repo: https://github.com/muneeb-dotcom/plp1-pmd-digital-twin

## 3. Therapeutic Design (Module B)
Mechanism-matched intervention map built across all 3 classes: AAV replacement (loss-of-function, dosage-limited MBP promoter to avoid overshoot), ASO knockdown (duplication, 82 candidate 3'UTR target sites), and chaperone/PERK-inhibition (misfolding, target prioritized via Open Targets, 6 compounds docked against the Gly246Trp/Leu31 pocket, best -7.7 kcal/mol). CNS-penetrance filtering found no compound fully brain-ready -- a real gap, not a dead end.

# Module C -- Final Synthesis (Phase C7)

**Scope note:** this covers Module C (Parts 5-7) only, assembled directly from the result files this module produced. Sections 1-3 of the full project manuscript (PLP1 mechanistic heterogeneity / Module A cell-state clustering / Module B intervention matching) come from those modules' own outputs and are not re-derived here -- merge them in from their result files before submission.

## 4. Restoration and gating (Parts 5-6)
Of 9 mechanism x intervention cells tested, 1 had real measured data (misfolding_confirmed x chaperone/PERK-inhibition, GSE277705).
Restoration metrics: 3 of 6 were computable from real data (survival_index, myelin_recovery, stress_resolution). maturation_recovery, plp1_normalisation, and ol_identity had no measurable equivalent in this bulk dataset.
Gated vs naive myelination under PERK inhibition: gated delta = 0.176, naive delta = 0.075 (both positive -- PERK inhibition rescues survival and myelination together, unlike the published Ro 25-6981 dissociation it was compared against in Test 4).

## 5. Held-out validation (Part 7)

| Test | Result | Status |
|---|---|---|
| 1. Mechanism classification | 0.562 vs 0.5 chance (n=12) | FAILED |
| 2. Severity correlation | rho=-0.034, p=0.92 (n=12) | FAILED |
| 3. Intervention response | no data exists | NOT MEASURABLE |
| 4. Dissociation | different real intervention, no match | NOT REPRODUCED |

0 of 4 tests passed. Reported as-is: a validation section with no failures is less credible than one with honestly characterized negative results.

Full test-by-test detail: `results/validation_report.md`.
Failure characterization, including the PLP1-null displacement/severity mismatch: `results/failure_analysis.md`.

## 6. Falsifiable experiments proposed
Three predictions, each with a stated falsification condition, in `results/experimental_plan.md`:
1. PLP1-null vs point-mutant oligodendrocytes fail via different mechanisms (ER stress vs not)
2. PLP1 copy-number response is non-monotonic (triplication showed lower displacement than duplication despite higher clinical severity)
3. PERK inhibition rescues myelination where Ro 25-6981 does not, in a head-to-head assay

## Headline result of Module C
0 of 4 pre-specified validation tests passed. The single most informative finding is not a pass -- it's that transcriptional displacement-from-healthy and clinical severity appear to track different underlying biology (PMD12, PLP1-null, clinically mildest but transcriptionally most displaced of all 12 lines). That result survived a circularity correction (PLP1 removed from its own scoring panel) and is proposed as Experimental Prediction 1.