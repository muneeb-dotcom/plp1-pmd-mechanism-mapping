import pandas as pd
import numpy as np
from sklearn.metrics import balanced_accuracy_score

state = pd.read_csv("results/nevin_state_vectors.csv")
meta = pd.read_csv("data/heldout/gse96049/nevin_line_metadata.csv")

state["line"] = state["sample"].str.replace(r"[ab]$", "", regex=True)
line_state = state.groupby("line")[["myelin_output","apoptosis","upr_signaling"]].mean()

pmd = line_state.loc[[l for l in line_state.index if l.startswith("PMD")]]
pmd = pmd.merge(meta, left_index=True, right_on="line").set_index("line")

pmd["mechanism_binary"] = pmd["mutation_type"].apply(
    lambda m: "point_mutation" if m == "point_mutation" else "structural_variant"
)

def loo_accuracy(pmd_subset):
    X = pmd_subset[["myelin_output","apoptosis","upr_signaling"]].values
    y = pmd_subset["mechanism_binary"].values
    baseline = pd.Series(y).value_counts(normalize=True).max()
    preds = []
    for i in range(len(y)):
        mask = np.ones(len(y), dtype=bool); mask[i] = False
        centroids = {c: X[mask][y[mask] == c].mean(axis=0) for c in set(y[mask])}
        dists = {c: np.linalg.norm(X[i] - cen) for c, cen in centroids.items()}
        preds.append(min(dists, key=dists.get))
    acc = balanced_accuracy_score(y, preds)
    return acc, baseline, preds

acc, baseline_acc, preds = loo_accuracy(pmd)
pmd["predicted_binary"] = preds
print(f"ORIGINAL (n={len(pmd)}): balanced accuracy={acc:.3f} vs baseline={baseline_acc:.3f}")

pmd_sens = pmd.drop("PMD9")
acc_s, baseline_s, _ = loo_accuracy(pmd_sens)
print(f"SENSITIVITY, PMD9 excluded (n={len(pmd_sens)}): balanced accuracy={acc_s:.3f} vs baseline={baseline_s:.3f}")

pmd[["mutation_type","mechanism_binary","predicted_binary"]].to_csv("results/nevin_test1_binary.csv")