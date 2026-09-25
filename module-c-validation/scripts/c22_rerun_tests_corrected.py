import pandas as pd
import numpy as np
import glob, os, re
from scipy.stats import spearmanr
from sklearn.metrics import balanced_accuracy_score

MODULES = {
    "myelin_output": ["MBP","MOG","MAG","CNP","MAL","UGT8","MOBP"],  # PLP1 removed
    "apoptosis": ["CASP3","BAX","BCL2L11","TP53","CDKN1A"],
    "upr_signaling": ["ATF4","DDIT3","EIF2AK3","XBP1","ATF6"],
}

files = glob.glob("data/heldout/gse96049/extracted/*.fpkm_tracking.gz")
rows = []
for f in files:
    m = re.search(r"GSM\d+_([A-Za-z0-9]+)\.fpkm_tracking", os.path.basename(f))
    sample = m.group(1)
    d = pd.read_csv(f, sep="\t")[["gene_short_name", "FPKM"]]
    d.columns = ["gene", "FPKM"]
    d = d.groupby("gene", as_index=False).mean()
    log_fpkm = np.log2(d.set_index("gene")["FPKM"] + 1)
    row = {"sample": sample}
    for module, genes in MODULES.items():
        present = [g for g in genes if g in log_fpkm.index]
        row[module] = log_fpkm.loc[present].mean()
    rows.append(row)

state = pd.DataFrame(rows)
state["line"] = state["sample"].str.replace(r"[ab]$", "", regex=True)
line_state = state.groupby("line")[["myelin_output","apoptosis","upr_signaling"]].mean()

meta = pd.read_csv("data/heldout/gse96049/nevin_line_metadata.csv")
pmd = line_state.loc[[l for l in line_state.index if l.startswith("PMD")]]
nc = line_state.loc[[l for l in line_state.index if l.startswith("NC")]]
pmd = pmd.merge(meta, left_index=True, right_on="line").set_index("line")

z = (line_state - line_state.mean()) / line_state.std()
nc_centroid_z = z.loc[nc.index].mean()
pmd_z = z.loc[pmd.index]
pmd["displacement"] = np.sqrt(((pmd_z - nc_centroid_z) ** 2).sum(axis=1))

severity_rank = {"mild": 1, "moderate": 2, "severe": 3}
pmd["severity_rank"] = pmd["severity"].map(severity_rank)
rho, p = spearmanr(pmd["displacement"], pmd["severity_rank"])
print("TEST 2 (corrected, PLP1 removed from panel)")
print(pmd[["mutation_type","severity","displacement"]].to_string())
print(f"Spearman rho={rho:.3f}, p={p:.4f}, n={len(pmd)}\n")

pmd["mechanism_binary"] = pmd["mutation_type"].apply(
    lambda m: "point_mutation" if m == "point_mutation" else "structural_variant"
)
X = pmd[["myelin_output","apoptosis","upr_signaling"]].values
y = pmd["mechanism_binary"].values
preds = []
for i in range(len(y)):
    mask = np.ones(len(y), dtype=bool); mask[i] = False
    centroids = {c: X[mask][y[mask] == c].mean(axis=0) for c in set(y[mask])}
    dists = {c: np.linalg.norm(X[i] - cen) for c, cen in centroids.items()}
    preds.append(min(dists, key=dists.get))
acc = balanced_accuracy_score(y, preds)
print("TEST 1 (corrected, PLP1 removed from panel)")
print(f"balanced accuracy={acc:.3f} vs chance=0.5, n={len(y)}")

pmd.to_csv("results/nevin_tests_corrected.csv")