import pandas as pd
import matplotlib.pyplot as plt

perk = pd.read_csv("results/gate_math_validation.csv", index_col=0)
rest = pd.read_csv("results/restoration_metrics.csv")
rest = rest[rest["status"] == "measured"]

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

x = ["active", "KO"]
width = 0.35
pos = range(len(x))
axes[0].bar([p - width/2 for p in pos], perk.loc[x, "naive_myelination"], width, label="naive", color="#999999")
axes[0].bar([p + width/2 for p in pos], perk.loc[x, "gated_myelination"], width, label="gated", color="#2a6f97")
axes[0].set_xticks(list(pos)); axes[0].set_xticklabels(x)
axes[0].set_ylabel("myelination score")
axes[0].set_title("Naive vs gated myelination\n(GSE277705, PERK KO)")
axes[0].legend()

axes[1].bar(rest["metric"], rest["raw_shift"], color="#2a6f97")
axes[1].axhline(0, color="black", linewidth=0.8)
axes[1].set_ylabel("KO - active (z-score)")
axes[1].set_title("Restoration metrics\n(measured only)")
axes[1].tick_params(axis="x", rotation=25)

plt.tight_layout()
plt.savefig("figures/c16_gate_and_restoration.png", dpi=200)
print("saved figures/c16_gate_and_restoration.png")