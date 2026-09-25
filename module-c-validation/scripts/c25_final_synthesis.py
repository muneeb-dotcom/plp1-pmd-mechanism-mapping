import pandas as pd

def read(path):
    with open(path) as f:
        return f.read()

restoration = pd.read_csv("results/restoration_metrics.csv")
restoration_measured = restoration[restoration["status"] == "measured"]
gate = pd.read_csv("results/gate_math_validation.csv", index_col=0)
predicted = pd.read_csv("results/predicted_responses.csv")
validation_report = read("results/validation_report.md")
failure_analysis = read("results/failure_analysis.md")
experimental_plan = read("results/experimental_plan.md")

lines = []
lines.append("# Module C -- Final Synthesis (Phase C7)\n")
lines.append("**Scope note:** this covers Module C (Parts 5-7) only, assembled directly from the result "
              "files this module produced. Sections 1-3 of the full project manuscript (PLP1 mechanistic "
              "heterogeneity / Module A cell-state clustering / Module B intervention matching) come from "
              "those modules' own outputs and are not re-derived here -- merge them in from their result "
              "files before submission.\n")

lines.append("## 4. Restoration and gating (Parts 5-6)")
lines.append(f"Of {len(predicted)} mechanism x intervention cells tested, "
              f"{(predicted['status']=='measured').sum()} had real measured data "
              f"(misfolding_confirmed x chaperone/PERK-inhibition, GSE277705).")
lines.append(f"Restoration metrics: {len(restoration_measured)} of {len(restoration)} were computable from "
              "real data (survival_index, myelin_recovery, stress_resolution). maturation_recovery, "
              "plp1_normalisation, and ol_identity had no measurable equivalent in this bulk dataset.")
lines.append(f"Gated vs naive myelination under PERK inhibition: gated delta = "
              f"{gate.loc['KO','gated_myelination'] - gate.loc['active','gated_myelination']:.3f}, "
              f"naive delta = {gate.loc['KO','naive_myelination'] - gate.loc['active','naive_myelination']:.3f} "
              "(both positive -- PERK inhibition rescues survival and myelination together, unlike the "
              "published Ro 25-6981 dissociation it was compared against in Test 4).\n")

lines.append("## 5. Held-out validation (Part 7)")
lines.append(validation_report.split("## Summary")[1] if "## Summary" in validation_report
             else "See validation_report.md")
lines.append("")
lines.append("Full test-by-test detail: `results/validation_report.md`.")
lines.append("Failure characterization, including the PLP1-null displacement/severity mismatch: "
              "`results/failure_analysis.md`.\n")

lines.append("## 6. Falsifiable experiments proposed")
lines.append("Three predictions, each with a stated falsification condition, in `results/experimental_plan.md`:")
lines.append("1. PLP1-null vs point-mutant oligodendrocytes fail via different mechanisms (ER stress vs not)")
lines.append("2. PLP1 copy-number response is non-monotonic (triplication showed lower displacement than "
              "duplication despite higher clinical severity)")
lines.append("3. PERK inhibition rescues myelination where Ro 25-6981 does not, in a head-to-head assay\n")

lines.append("## Headline result of Module C")
lines.append("0 of 4 pre-specified validation tests passed. The single most informative finding is not a "
              "pass -- it's that transcriptional displacement-from-healthy and clinical severity appear to "
              "track different underlying biology (PMD12, PLP1-null, clinically mildest but transcriptionally "
              "most displaced of all 12 lines). That result survived a circularity correction (PLP1 removed "
              "from its own scoring panel) and is proposed as Experimental Prediction 1.")

with open("results/final_synthesis.md", "w") as f:
    f.write("\n".join(lines))

print("\n".join(lines))