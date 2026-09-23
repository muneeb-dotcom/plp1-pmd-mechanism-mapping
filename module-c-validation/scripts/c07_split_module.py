import pandas as pd
import numpy as np
import json

with open("data/gene_symbol_to_ensembl.json") as f:
    sym2ens = json.load(f)

files = {
    "neg1": "GSM8527875_Jimpy_perk_Flox_Cre_negative_1_S10.txt.gz",
    "neg2": "GSM8527876_Jimpy_perk_Flox_Cre_negative_2_S11.txt.gz",
    "neg3": "GSM8527877_Jimpy_perk_Flox_Cre_negative_3_S12.txt.gz",
    "neg4": "GSM8527878_Jimpy_perk_Flox_Cre_negative_4_S13.txt.gz",
    "pos1": "GSM8527879_Jimpy_perk_Flox_Cre_positive_1_S14.txt.gz",
    "pos2": "GSM8527880_Jimpy_perk_Flox_Cre_positive_2_S15.txt.gz",
    "pos3": "GSM8527881_Jimpy_perk_Flox_Cre_positive_3_S16.txt.gz",
    "pos4": "GSM8527882_Jimpy_perk_Flox_Cre_positive_4_S17.txt.gz",
}
base = "data/heldout/gse277705/extracted/"

counts = {}
for name, fname in files.items():
    df = pd.read_csv(base + fname, sep="\t", compression="gzip", index_col="Gene")
    counts[name] = df["Counts"]
mat = pd.DataFrame(counts)
cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

neg_cols = [c for c in log_cpm.columns if c.startswith("neg")]
pos_cols = [c for c in log_cpm.columns if c.startswith("pos")]

SPLIT_MODULES = {
    "upr_signaling": ["Atf4", "Ddit3", "Eif2ak3", "Xbp1", "Atf6"],
    "chaperone_capacity": ["Hspa5", "Dnajb9", "Pdia4", "Pdia6", "Manf", "Herpud1", "Hyou1"],
}

results = []
for module, genes in SPLIT_MODULES.items():
    ens_ids = [sym2ens[g] for g in genes if g in sym2ens]
    present = [e for e in ens_ids if e in log_cpm.index]
    neg_mean = log_cpm.loc[present, neg_cols].mean(axis=0).mean()
    pos_mean = log_cpm.loc[present, pos_cols].mean(axis=0).mean()
    results.append({
        "module": module, "genes_found": f"{len(present)}/{len(genes)}",
        "jimpy_PERK_active": neg_mean, "jimpy_PERK_knockout": pos_mean,
        "delta": pos_mean - neg_mean,
    })

df = pd.DataFrame(results)
print(df.to_string(index=False))
df.to_csv("results/split_module_validation.csv", index=False)