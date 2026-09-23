import pandas as pd
import numpy as np
import json
import glob

with open("data/gene_symbol_to_ensembl.json") as f:
    sym2ens = json.load(f)
with open("../plp1-twin/results/state_modules.json") as f:
    ref = json.load(f)
MODULE_GENES = ref["genes"]

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
print(f"Count matrix: {mat.shape}")

cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

neg_cols = [c for c in log_cpm.columns if c.startswith("neg")]
pos_cols = [c for c in log_cpm.columns if c.startswith("pos")]

results = []
for module, genes in MODULE_GENES.items():
    ens_ids = [sym2ens[g] for g in genes if g in sym2ens]
    present = [e for e in ens_ids if e in log_cpm.index]
    if not present:
        continue
    neg_mean = log_cpm.loc[present, neg_cols].mean(axis=0).mean()
    pos_mean = log_cpm.loc[present, pos_cols].mean(axis=0).mean()
    results.append({
        "module": module,
        "genes_found": f"{len(present)}/{len(genes)}",
        "jimpy_PERK_active_mean": neg_mean,
        "jimpy_PERK_knockout_mean": pos_mean,
        "delta_KO_minus_active": pos_mean - neg_mean,
    })

df_res = pd.DataFrame(results)
print(df_res.to_string(index=False))
df_res.to_csv("results/perk_knockout_validation.csv", index=False)

print("\nPRE-REGISTERED PREDICTION CHECK:")
er = df_res[df_res["module"] == "er_stress"]["delta_KO_minus_active"].values
myelin = df_res[df_res["module"] == "myelin_output"]["delta_KO_minus_active"].values
print(f"  er_stress delta (expect NEGATIVE): {er}")
print(f"  myelin_output delta (expect POSITIVE): {myelin}")