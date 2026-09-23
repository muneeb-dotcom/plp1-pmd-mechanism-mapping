import requests
import pandas as pd
import time

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

TARGETS = {
    "PLP1": "ENSG00000123560",
    "EIF2AK3": "ENSG00000172071",
    "ATF4": "ENSG00000128272",
    "DDIT3": "ENSG00000175197",
    "HSPA5": "ENSG00000044574",
    "ERN1": "ENSG00000178607",
    "ATF6": "ENSG00000118217",
}

QUERY = """
query($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    approvedSymbol
    approvedName
    isEssential
    tractability { label modality value }
    geneticConstraint { constraintType oe oeUpper upperBin6 }
    safetyLiabilities { event datasource }
  }
}
"""

def summarize_tractability(tract_list):
    sm = {t["label"]: t["value"] for t in tract_list if t["modality"] == "SM"}
    ab = {t["label"]: t["value"] for t in tract_list if t["modality"] == "AB"}
    return {
        "SM_druggable_family": sm.get("Druggable Family", False),
        "SM_has_ligand":       sm.get("High-Quality Ligand", False),
        "SM_has_pocket":       sm.get("High-Quality Pocket", False) or sm.get("Med-Quality Pocket", False),
        "SM_approved_drug":    sm.get("Approved Drug", False),
        "AB_feasible":         ab.get("UniProt loc high conf", False),
    }

rows = []
for symbol, ensembl_id in TARGETS.items():
    r = requests.post(API_URL, json={"query": QUERY, "variables": {"ensemblId": ensembl_id}}, timeout=30)
    r.raise_for_status()
    data = r.json()["data"]["target"]

    tract = summarize_tractability(data["tractability"])
    lof_constraint = next((c for c in data["geneticConstraint"] if c["constraintType"] == "lof"), {})
    safety_events = [s["event"] for s in data["safetyLiabilities"]]

    rows.append({
        "symbol": data["approvedSymbol"],
        "name": data["approvedName"],
        "is_essential": data["isEssential"],
        **tract,
        "lof_oe": lof_constraint.get("oe"),
        "lof_constraint_bin": lof_constraint.get("upperBin6"),
        "safety_events": "; ".join(safety_events) if safety_events else "none reported",
    })
    time.sleep(0.3)  # be polite to the public API

df = pd.DataFrame(rows)
df["tractable_smallmolecule"] = df["SM_druggable_family"] & (df["SM_has_ligand"] | df["SM_has_pocket"])
df.to_csv("../results/target_priority.csv", index=False)

print(df[["symbol", "is_essential", "tractable_smallmolecule", "lof_oe", "safety_events"]].to_string(index=False))