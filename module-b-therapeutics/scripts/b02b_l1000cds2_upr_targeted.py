import pandas as pd
import requests
import json

de = pd.read_csv("results/plp1_jimpy_deseq_results.csv", index_col=0)

UPR_PANEL = ["Atf4", "Atf6", "Ddit3", "Xbp1", "Hspa5", "Herpud1", "Eif2ak3",
             "Ern1", "Dnajb9", "Edem1", "Pdia4", "Pdia6", "Hyou1", "Sec61a1",
             "Calr", "Canx", "Manf", "Sel1l", "Derl1"]

present = de.reindex(UPR_PANEL).dropna(subset=["t"])
print(f"{len(present)}/{len(UPR_PANEL)} UPR panel genes found with t-statistics")

up_in_disease   = present[present["t"] > 0].index.tolist()
down_in_disease = present[present["t"] < 0].index.tolist()
print(f"Up in disease (jimpy):   {up_in_disease}")
print(f"Down in disease (jimpy): {down_in_disease}")

# Reversal query: a correcting compound should push the disease-up genes
# DOWN and the disease-down genes UP -- so disease-up genes go in dnGenes
# (we want them decreased) and disease-down genes go in upGenes.
payload = {
    "data": {
        "upGenes": [g.upper() for g in down_in_disease],
        "dnGenes": [g.upper() for g in up_in_disease],
    },
    "config": {
        "aggravate": False,
        "searchMethod": "geneSet",
        "share": False,
        "combination": False,
        "db-version": "latest",
    },
}

print("\nSubmitting UPR-targeted reversal query to L1000CDS2...")
r = requests.post(
    "https://maayanlab.cloud/L1000CDS2/query",
    data=json.dumps(payload),
    headers={"content-type": "application/json"},
    timeout=60,
)
r.raise_for_status()
result = r.json()

if "err" in result:
    print("API error:", result["err"])
else:
    hits = pd.DataFrame(result["topMeta"])
    cols = ["pert_desc", "pert_id", "cell_id", "pert_dose", "pert_dose_unit", "score"]
    hits = hits[[c for c in cols if c in hits.columns]]
    hits.to_csv("results/l1000cds2_upr_targeted_hits.csv", index=False)

    print(f"\n{len(hits)} UPR-targeted reversal candidates returned")
    print(hits.head(20).to_string(index=False))

    targets_of_interest = ["ISRIB", "SALUBRINAL", "GSK2606414", "TUNICAMYCIN",
                           "4-PBA", "PHENYLBUTYRATE", "TAUROURSODEOXYCHOLIC"]
    print("\nChecking for known UPR-modulating compounds:")
    for t in targets_of_interest:
        match = hits[hits["pert_desc"].astype(str).str.upper().str.contains(t, na=False)]
        if len(match):
            print(f"  FOUND: {t}")
            print(match.to_string(index=False))
        else:
            print(f"  not found: {t}")