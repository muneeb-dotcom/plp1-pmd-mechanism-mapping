# Module C -- Held-Out Validation Report

Four pre-specified tests (Module C, Phase C4). Results reported as run, including failures.

## Test 1 -- Mechanism classification (Nevin et al. 2017, GSE96049)
**Criterion:** balanced accuracy vs chance (0.5 for binary).
**Result:** balanced accuracy = 0.562 vs chance = 0.5, n = 12. **FAILED** -- margin of 0.062 on n=12 is not distinguishable from noise.
An uncorrected first pass scored 0.312 (worse than chance); that panel included PLP1 itself as a marker gene, which is circular for the PLP1-deletion line (PMD12, PLP1_FPKM = 0 by construction). Corrected panel excludes PLP1.
Binary split (point_mutation n=8 vs structural_variant n=4) used after the original 5-class split produced singleton classes that are structurally unclassifiable under leave-one-out.
**Interpretation:** the state vector does not separate point mutations from structural variants in this cohort, before or after correction. Real negative result.

## Test 2 -- Severity correlation (Nevin et al. 2017, GSE96049)
**Criterion:** Spearman rho between predicted displacement and clinical severity rank.
**Result:** rho = -0.034, p = 0.9165, n = 12. **FAILED** (no correlation).
**Correction note:** myelin_output originally included PLP1 as a marker gene -- circular for PMD12 (full deletion). Removing PLP1 reduced PMD12's displacement from 4.54 to 2.88, but it remains the single largest displacement of all 12 lines despite being the clinically mildest. Real, non-circular finding: PLP1-null is transcriptionally distinct from wild-type in this state space even though loss-of-function is clinically milder than toxic gain-of-function point mutations -- displacement and clinical severity may track different biology. Flagged for C5.

## Test 3 -- Intervention response (Elitt et al. 2020, Nature, PMID 32610343)
**Criterion:** directional concordance between predicted and reported response to PLP1-ASO.
**Result: NOT_MEASURABLE.** No RNA-seq deposited for this study (readouts are histology, motor/respiratory function, lifespan). No fitted twin in this pipeline produces a scalar prediction to score against those readouts.
Reported direction only: oligodendrocyte numbers up, myelination up, motor/respiratory improved, lifespan extended -- survival AND myelination co-improve.

## Test 4 -- Survival/myelination dissociation (Elitt et al. 2018, PMID 30146490)
**Criterion:** binary -- does a real tested intervention reproduce survival-rescued-without-myelination?
**Result:** PERK inhibition rescues BOTH survival and myelination -- does not reproduce the Ro 25-6981 dissociation. Different mechanism, different real outcome.

## Summary
| Test | Result | Status |
|---|---|---|
| 1. Mechanism classification | 0.562 vs 0.5 chance (n=12) | FAILED |
| 2. Severity correlation | rho=-0.034, p=0.92 (n=12) | FAILED |
| 3. Intervention response | no data exists | NOT MEASURABLE |
| 4. Dissociation | different real intervention, no match | NOT REPRODUCED |

0 of 4 tests passed. Reported as-is: a validation section with no failures is less credible than one with honestly characterized negative results.