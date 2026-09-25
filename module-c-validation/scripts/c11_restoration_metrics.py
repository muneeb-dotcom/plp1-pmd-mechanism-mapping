import pandas as pd

perk = pd.read_csv("results/gate_math_validation.csv", index_col=0)

METRICS = {
    "survival_index":     "apoptosis",
    "myelin_recovery":    "myelin_output",
    "stress_resolution":  "upr_signaling",
}

rows = []
for label, col in METRICS.items():
    active = perk.loc["active", col]
    ko = perk.loc["KO", col]
    rows.append({"metric": label, "module": col, "status": "measured",
                  "active": active, "KO": ko, "raw_shift": ko - active})

for extra, why in [
    ("maturation_recovery", "no maturation axis exists for bulk whole-tissue samples"),
    ("plp1_normalisation", "no PLP1 dosage measurement in this dataset"),
    ("ol_identity", "no single-cell centroid data in this dataset"),
]:
    rows.append({"metric": extra, "module": None, "status": f"no_measured_data: {why}",
                  "active": None, "KO": None, "raw_shift": None})

df = pd.DataFrame(rows)
df.to_csv("results/restoration_metrics.csv", index=False)
print(df.to_string(index=False))