import pandas as pd
import numpy as np
import json

with open("data/gene_symbol_to_ensembl.json") as f:
    sym2ens = json.load(f)

files = {
    "neg1": "GSM8527867_Perk_Group_Cre_negative_1_S22.txt.gz",
    "neg2": "GSM8527868_Perk_Group_Cre_negative_2_S23.txt.gz",
    "neg3": "GSM8527869_Perk_Group_Cre_negative_3_S24.txt.gz",
    "neg4": "GSM8527870_Perk_Group_Cre_negative_4_S25.txt.gz",
    "pos1": "GSM8527871_Perk_Group_Cre_positive_1_S18.txt.gz",
    "pos2": "GSM8527872_Perk_Group_Cre_positive_2_S19.txt.gz",
    "pos3": "GSM8527873_Perk_Group_Cre_positive_3_S20.txt.gz",
    "pos4": "GSM8527874_Perk_Group_Cre_positive_4_S21.txt.gz",
}
base = "data/heldout/gse277705/extracted/"

counts = {n: pd.read_csv(base + f, sep="\t", compression="gzip", index_col="Gene")["Counts"]
          for n, f in files.items()}
mat = pd.DataFrame(counts)
cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
log_cpm = np.log2(cpm + 1)

MODULES = {
    "myelin_output": ["Mbp","Mog","Mag","Cnp","Mal","Plp1","Ugt8a","Mobp"],
    "apoptosis": ["Casp3","Bax","Bcl2l11","Trp53","Cdkn1a"],
    "upr_signaling": ["Atf4","Ddit3","Eif2ak3","Xbp1","Atf6"],
}

neg_cols = [c for c in log_cpm.columns if c.startswith("neg")]
pos_cols = [c for c in log_cpm.columns if c.startswith("pos")]

rows = []
for module, genes in MODULES.items():
    ens_ids = [sym2ens[g] for g in genes if g in sym2ens]
    present = [e for e in ens_ids if e in log_cpm.index]
    neg_mean = log_cpm.loc[present, neg_cols].mean(axis=0).mean()
    pos_mean = log_cpm.loc[present, pos_cols].mean(axis=0).mean()
    rows.append({"module": module, "background_active": neg_mean,
                 "background_PERK_KO": pos_mean, "delta_background_only": pos_mean - neg_mean})

df = pd.DataFrame(rows)
print("PERK-KO effect in NON-jimpy (wild-type genetic background):")
print(df.round(3).to_string(index=False))

jimpy_arm = pd.read_csv("results/perk_knockout_validation.csv") if False else None
print("\nCompare against jimpy-background PERK-KO delta (from c04/c07):")
print("  myelin_output delta was +1.005 in jimpy background")
print("  upr_signaling delta was -0.039 in jimpy background")
print("\nIf background-only deltas are much smaller, the myelin rescue is jimpy/disease-specific.")
print("If similar in magnitude, PERK-KO boosts myelination generally, independent of disease.")

df.to_csv("results/perk_background_specificity_check.csv", index=False)