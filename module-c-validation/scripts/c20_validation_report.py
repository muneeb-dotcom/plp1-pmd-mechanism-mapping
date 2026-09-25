import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import balanced_accuracy_score
import numpy as np

nevin = pd.read_csv("results/nevin_tests_corrected.csv")
t3 = pd.read_csv("results/test3_intervention_response.csv")
t4 = pd.read_csv("results/myelination_predictions.csv")

# recompute stats directly from the corrected data rather than hardcoding
rho, p = spearmanr(nevin["displacement"], nevin["severity_rank"])

X = nevin[["myelin_output","apoptosis","upr_signaling"]].values
y = nevin["mechanism_binary"].values
preds = []
for i in range(len(y)):
    mask = np.ones(len(y), dtype=bool); mask[i] = False
    centroids = {c: X[mask][y[mask] == c].mean(axis=0) for c in set(y[mask])}
    dists = {c: np.linalg.norm(X[i] - cen) for c, cen in centroids.items()}
    preds.append(min(dists, key=dists.get))
acc = balanced_accuracy_score(y, preds)

lines = []
lines.append("# Module C -- Held-Out Validation Report\n")
lines.append("Four pre-specified tests (Module C, Phase C4). Results reported as run, including failures.\n")

lines.append("## Test 1 -- Mechanism classification (Nevin et al. 2017, GSE96049)")
lines.append("**Criterion:** balanced accuracy vs chance (0.5 for binary).")
lines.append(f"**Result:** balanced accuracy = {acc:.3f} vs chance = 0.5, n = {len(y)}. "
              f"**FAILED** -- margin of {acc-0.5:.3f} on n={len(y)} is not distinguishable from noise.")
lines.append("An uncorrected first pass scored 0.312 (worse than chance); that panel included PLP1 itself as "
              "a marker gene, which is circular for the PLP1-deletion line (PMD12, PLP1_FPKM = 0 by "
              "construction). Corrected panel excludes PLP1.")
lines.append("Binary split (point_mutation n=8 vs structural_variant n=4) used after the original 5-class "
              "split produced singleton classes that are structurally unclassifiable under leave-one-out.")
lines.append("**Interpretation:** the state vector does not separate point mutations from structural "
              "variants in this cohort, before or after correction. Real negative result.\n")

lines.append("## Test 2 -- Severity correlation (Nevin et al. 2017, GSE96049)")
lines.append("**Criterion:** Spearman rho between predicted displacement and clinical severity rank.")
lines.append(f"**Result:** rho = {rho:.3f}, p = {p:.4f}, n = {len(nevin)}. **FAILED** (no correlation).")
lines.append("**Correction note:** myelin_output originally included PLP1 as a marker gene -- circular for "
              "PMD12 (full deletion). Removing PLP1 reduced PMD12's displacement from 4.54 to 2.88, but it "
              "remains the single largest displacement of all 12 lines despite being the clinically mildest. "
              "Real, non-circular finding: PLP1-null is transcriptionally distinct from wild-type in this "
              "state space even though loss-of-function is clinically milder than toxic gain-of-function "
              "point mutations -- displacement and clinical severity may track different biology. "
              "Flagged for C5.\n")

lines.append("## Test 3 -- Intervention response (Elitt et al. 2020, Nature, PMID 32610343)")
lines.append("**Criterion:** directional concordance between predicted and reported response to PLP1-ASO.")
row = t3.iloc[0]
lines.append(f"**Result: {row['status'].upper()}.** No RNA-seq deposited for this study (readouts are "
              "histology, motor/respiratory function, lifespan). No fitted twin in this pipeline produces a "
              "scalar prediction to score against those readouts.")
lines.append(f"Reported direction only: {row['reported_direction']}.\n")

lines.append("## Test 4 -- Survival/myelination dissociation (Elitt et al. 2018, PMID 30146490)")
lines.append("**Criterion:** binary -- does a real tested intervention reproduce survival-rescued-without-myelination?")
row4 = t4.iloc[0]
lines.append(f"**Result:** {row4['result']}\n")

lines.append("## Summary")
lines.append("| Test | Result | Status |")
lines.append("|---|---|---|")
lines.append(f"| 1. Mechanism classification | {acc:.3f} vs 0.5 chance (n={len(y)}) | FAILED |")
lines.append(f"| 2. Severity correlation | rho={rho:.3f}, p={p:.2f} (n={len(nevin)}) | FAILED |")
lines.append("| 3. Intervention response | no data exists | NOT MEASURABLE |")
lines.append("| 4. Dissociation | different real intervention, no match | NOT REPRODUCED |")
lines.append("\n0 of 4 tests passed. Reported as-is: a validation section with no failures is less credible "
              "than one with honestly characterized negative results.")

with open("results/validation_report.md", "w") as f:
    f.write("\n".join(lines))

print("\n".join(lines))