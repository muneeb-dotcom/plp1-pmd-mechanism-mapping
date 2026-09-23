import pandas as pd
import requests
import json

de = pd.read_csv("results/plp1_jimpy_deseq_results.csv", index_col=0)
de = de.dropna(subset=["t"])

N = 150
up_genes   = de.nlargest(N, "t").index.tolist()
down_genes = de.nsmallest(N, "t").index.tolist()

# Mouse -> human symbol approximation: simple uppercase conversion.
# This is a first-pass approximation, not true ortholog mapping (same
# caveat noted for the digital twin's species bridge) -- most mouse
# gene symbols do uppercase directly to their human HGNC equivalent,
# but gene-family members with different naming conventions can be missed.
up_human   = [g.upper() for g in up_genes]
down_human = [g.upper() for g in down_genes]

payload = {
    "data": {
        "upGenes": up_human,
        "dnGenes": down_human,
    },
    "config": {
        "aggravate": False,   # False = search for REVERSAL compounds (what we want)
        "searchMethod": "geneSet",
        "share": False,
        "combination": False,
        "db-version": "latest",
    },
}

print("Submitting query to L1000CDS2 (reversal mode)...")
r = requests.post(
    "https://maayanlab.cloud/L1000CDS2/query",
    data=json.dumps(payload),
    headers={"content-type": "application/json"},
    timeout=60,
)
r.raise_for_status()
result = r.json()

if "err" in result:
    print("API returned an error:", result["err"])
else:
    hits = pd.DataFrame(result["topMeta"])
    cols = ["pert_desc", "pert_id", "cell_id", "pert_dose", "pert_dose_unit", "score"]
    hits = hits[[c for c in cols if c in hits.columns]]
    hits.to_csv("results/l1000cds2_reversal_hits.csv", index=False)

    print(f"\n{len(hits)} reversal candidates returned")
    print(hits.head(20).to_string(index=False))

    targets_of_interest = ["ISRIB", "SALUBRINAL", "GSK2606414", "TUNICAMYCIN"]
    print("\nChecking for known UPR-modulating compounds in the results:")
    for t in targets_of_interest:
        match = hits[hits["pert_desc"].astype(str).str.upper().str.contains(t, na=False)]
        if len(match):
            print(f"  FOUND: {t}")
            print(match.to_string(index=False))
        else:
            print(f"  not found: {t}")