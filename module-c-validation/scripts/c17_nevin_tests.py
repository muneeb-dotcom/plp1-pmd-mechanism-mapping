import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import balanced_accuracy_score

state = pd.read_csv("results/nevin_state_vectors.csv")
meta = pd.read_csv("data/heldout/gse96049/nevin_line_metadata.csv")

state["line"] = state["sample"].str.replace(r"[ab]$", "", regex=True)
line_state = state.groupby("line")[["myelin_output","apoptosis","upr_signaling"]].mean()

pmd_lines = [l for l in line_state.index if l.startswith("PMD")]
nc_lines = [l for l in line_state.index if l.startswith("NC")]

z = (line_state - line_state.mean()) / line_state.std()
nc_centroid_z = z.loc[nc_lines].mean()

pmd = line_state.loc[pmd_lines].merge(meta, left_index=True, right_on="line").set_index("line")
pmd_z = z.loc[pmd_lines]
pmd["displacement"] = np.sqrt(((pmd_z - nc_centroid_z) ** 2).sum(axis=1))

severity_rank = {"mild": 1, "moderate": 2, "severe": 3}
pmd["severity_rank"] = pmd["severity"].map(severity_rank)
rho, p = spearmanr(pmd["displacement"], pmd["severity_rank"])
print("TEST 2 -- severity correlation")
print(pmd[["mutation_type","severity","displacement"]].to_string())
print(f"Spearman rho={rho:.3f}, p={p:.4f}, n={len(pmd)}\n")

X = pmd[["myelin_output","apoptosis","upr_signaling"]].values
y = pmd["mutation_type"].values
baseline_acc = pd.Series(y).value_counts(normalize=True).max()

preds = []
for i in range(len(y)):
    mask = np.ones(len(y), dtype=bool); mask[i] = False
    centroids = {c: X[mask][y[mask] == c].mean(axis=0) for c in set(y[mask])}
    dists = {c: np.linalg.norm(X[i] - cen) for c, cen in centroids.items()}
    preds.append(min(dists, key=dists.get))

acc = balanced_accuracy_score(y, preds)
print("TEST 1 -- mechanism classification (leave-one-out nearest centroid)")
print(f"balanced accuracy={acc:.3f} vs majority baseline={baseline_acc:.3f}, n={len(y)}")
print("NOTE: point_mutation n=8, all other classes n=1 -- classification is unreliable at this sample size")

pmd["predicted_mechanism"] = preds
pmd[["mutation_type","severity","displacement","predicted_mechanism"]].to_csv("results/nevin_validation_tests.csv")