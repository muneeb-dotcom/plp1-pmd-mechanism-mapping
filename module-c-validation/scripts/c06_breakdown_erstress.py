import pandas as pd
import numpy as np
import json

with open("data/gene_symbol_to_ensembl.json") as f:
    sym2ens = json.load(f)
with open("../plp1-twin/results/state_modules.json") as f:
    ref = json.load(f)

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

# rough UPR branch classification
perk_branch = ["Atf4", "Ddit3", "Eif2ak3"]
atf6_branch = ["Atf6", "Herpud1", "Hyou1"]
ire1_branch = ["Xbp1", "Ern1", "Edem1", "Sel1l", "Derl1"]
chaperones  = ["Hspa5", "Dnajb9", "Pdia4", "Pdia6", "Calr", "Canx", "Manf"]

for label, genes in [("PERK branch", perk_branch), ("ATF6 branch", atf6_branch),
                     ("IRE1 branch", ire1_branch), ("Chaperones", chaperones)]:
    print(f"\n{label}:")
    for g in genes:
        ens = sym2ens.get(g)
        if ens and ens in log_cpm.index:
            neg_m = log_cpm.loc[ens, neg_cols].mean()
            pos_m = log_cpm.loc[ens, pos_cols].mean()
            print(f"  {g}: active={neg_m:.3f}, KO={pos_m:.3f}, delta={pos_m-neg_m:+.3f}")
        else:
            print(f"  {g}: not in matrix")