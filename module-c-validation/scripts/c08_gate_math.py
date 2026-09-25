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

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

per_sample_scores = {}
for module, genes in MODULES.items():
    ens_ids = [sym2ens[g] for g in genes if g in sym2ens]
    present = [e for e in ens_ids if e in log_cpm.index]
    per_sample_scores[module] = log_cpm.loc[present].mean(axis=0)

scores_df = pd.DataFrame(per_sample_scores)

# within-dataset z-score (same platform, same study -- valid standardization)
z = (scores_df - scores_df.mean()) / scores_df.std()

# No maturation axis exists for bulk whole-tissue samples (single-cell-only
# concept from the atlas), so the survival term uses apoptosis alone here.
z["naive_myelination"] = sigmoid(-z["apoptosis"])
z["gated_myelination"] = (
    sigmoid(-z["apoptosis"]) * sigmoid(z["myelin_output"]) * sigmoid(-z["upr_signaling"])
)

z["group"] = ["active"]*4 + ["KO"]*4
summary = z.groupby("group")[["apoptosis","myelin_output","upr_signaling","naive_myelination","gated_myelination"]].mean()
print(z[["apoptosis","myelin_output","upr_signaling","naive_myelination","gated_myelination","group"]]
      .round(3).to_string())
print("\nGroup means:")
print(summary)
print("\nDelta (KO - active):")
print(summary.loc["KO"] - summary.loc["active"])
print("\nPaper says myelination increases under PERK KO.")
print("PREDICTION: gated_myelination delta should be POSITIVE and larger than naive delta.")

summary.to_csv("results/gate_math_validation.csv")