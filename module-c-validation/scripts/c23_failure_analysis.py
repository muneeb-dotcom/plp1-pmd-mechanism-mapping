import pandas as pd
import numpy as np

nevin = pd.read_csv("results/nevin_tests_corrected.csv")
t1 = pd.read_csv("results/nevin_test1_binary.csv")

lines = []
lines.append("# Module C -- Failure Analysis (Phase C5)\n")

lines.append("## 1. PLP1-null (PMD12) is the largest outlier despite being clinically mildest")
lines.append("Full PLP1 deletion (PMD12) shows the largest transcriptional displacement from wild-type "
              "of all 12 patient lines (2.88, vs range 0.81-2.30 for the rest), even with PLP1 itself "
              "excluded from the scoring panel to remove circularity. Clinically this line is the mildest "
              "(consistent with the literature: loss-of-function PLP1 mutations are generally milder than "
              "toxic gain-of-function point mutations). This means displacement-from-healthy, as scored here, "
              "is not a proxy for clinical severity -- it may instead be picking up a distinct axis (e.g. "
              "absence of a structural protein vs. presence of a misfolded one triggering ER stress). The two "
              "failure modes look transcriptionally similar in this state space but are clinically opposite.\n")

lines.append("## 2. Severity-rank vs displacement-rank mismatches")
nevin_sorted = nevin[["mutation_type","severity","severity_rank","displacement"]].copy()
nevin_sorted["displacement_rank"] = nevin_sorted["displacement"].rank()
nevin_sorted["rank_gap"] = (nevin_sorted["displacement_rank"] - nevin_sorted["severity_rank"] * (12/3)).abs()
worst = nevin_sorted.sort_values("rank_gap", ascending=False)
lines.append("Lines ranked by mismatch between clinical severity rank and displacement rank (largest gap first):")
lines.append(worst.to_string())
lines.append("")

lines.append("## 3. Mechanism misclassification pattern (Test 1, binary)")
confusion = t1.groupby(["mechanism_binary","predicted_binary"]).size().reset_index(name="n")
lines.append(confusion.to_string(index=False))
lines.append("All 4 structural_variant lines (duplication, triplication, full_deletion, partial_deletion) "
              "were misclassified as point_mutation. The state vector has no signal separating these classes "
              "-- structural variants are transcriptionally closer to the point-mutation centroid than to "
              "each other in this cohort.\n")

lines.append("## 4. Duplication vs triplication -- dosage non-linearity check")
dup = nevin.set_index("line").loc["PMD10"]
trip = nevin.set_index("line").loc["PMD11"]
lines.append(f"PMD10 (duplication, 2 copies): displacement = {dup['displacement']:.3f}, severity = {dup['severity']}")
lines.append(f"PMD11 (triplication, 3 copies): displacement = {trip['displacement']:.3f}, severity = {trip['severity']}")
ratio = trip['displacement'] / dup['displacement']
lines.append(f"Displacement ratio (3-copy / 2-copy) = {ratio:.2f}. If dosage effects were linear, tripling "
              "copy number relative to duplication (a 1.5x dosage step) would predict a proportionally larger "
              "displacement; a ratio far from ~1.5 indicates the model is not capturing dosage response "
              "linearly, but n=1 per class makes this a single-pair observation, not a trend.\n")

lines.append("## Summary")
lines.append("Two of four Part 7 tests failed outright (Test 1, Test 2); the most informative finding is not "
              "in the pass/fail table but in Section 1 above: the state vector conflates a mild loss-of-function "
              "mechanism with high transcriptional displacement, suggesting displacement-from-healthy and "
              "clinical severity are answering different biological questions and should not have been assumed "
              "to track each other going in.")

with open("results/failure_analysis.md", "w") as f:
    f.write("\n".join(lines))

print("\n".join(lines))