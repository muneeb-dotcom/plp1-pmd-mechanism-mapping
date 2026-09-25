import pandas as pd
import numpy as np

# The full 3x3 mechanism x intervention matrix. Only cells with a
# REAL measured dataset get a delta -- everything else stays honestly
# marked "no_measured_data" rather than being estimated or guessed.

perk_ko = pd.read_csv("results/gate_math_validation.csv", index_col=0)
naive_delta = perk_ko.loc["KO", "naive_myelination"] - perk_ko.loc["active", "naive_myelination"]
gated_delta = perk_ko.loc["KO", "gated_myelination"] - perk_ko.loc["active", "gated_myelination"]

MECHANISMS = ["loss_of_function", "duplication", "misfolding_confirmed"]
INTERVENTIONS = ["aav_replacement", "aso_knockdown", "chaperone_PERK_inhibition"]

rows = []
for mech in MECHANISMS:
    for iv in INTERVENTIONS:
        matched = (
            (mech == "loss_of_function" and iv == "aav_replacement") or
            (mech == "duplication" and iv == "aso_knockdown") or
            (mech == "misfolding_confirmed" and iv == "chaperone_PERK_inhibition")
        )
        if mech == "misfolding_confirmed" and iv == "chaperone_PERK_inhibition":
            status = "measured"
            naive_myelination_delta = naive_delta
            gated_myelination_delta = gated_delta
            source = "GSE277705 (real PERK knockout in jimpy mice)"
        else:
            status = "no_measured_data"
            naive_myelination_delta = None
            gated_myelination_delta = None
            source = None
        rows.append({
            "mechanism": mech, "intervention": iv, "matched": matched,
            "status": status,
            "naive_myelination_delta": naive_myelination_delta,
            "gated_myelination_delta": gated_myelination_delta,
            "source": source,
        })

df = pd.DataFrame(rows)
df.to_csv("results/predicted_responses.csv", index=False)
print("Phase C1 -- full intervention x mechanism matrix (honest coverage):\n")
print(df.to_string(index=False))

measured_count = (df["status"] == "measured").sum()
print(f"\n{measured_count} of {len(df)} cells have real measured data.")
print("This IS the finding: only the misfolding/PERK-inhibition cell is validated.")
print("All other cells rest on Module B literature precedent, not tested data.")