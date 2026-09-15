from Bio.PDB import PDBParser
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure("PLP1", "PLP1_alphafold_Repair.pdb")

ca_coords = {}
for residue in structure[0]["A"]:
    if "CA" in residue:
        ca_coords[residue.id[1]] = residue["CA"].coord

df = pd.read_csv("../../data/structures/plp1_structural_disruption_scores.csv")

xs, ys, zs, colors, labels = [], [], [], [], []
for _, row in df.iterrows():
    pos = int(row["protein_pos"])
    if pos in ca_coords:
        x, y, z = ca_coords[pos]
        xs.append(x); ys.append(y); zs.append(z)
        colors.append(row["ddG_kcal_mol"])
        labels.append(f"{row['ref_aa']}{pos}{row['alt_aa']}")

fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection="3d")

all_x = [c[0] for c in ca_coords.values()]
all_y = [c[1] for c in ca_coords.values()]
all_z = [c[2] for c in ca_coords.values()]
ax.plot(all_x, all_y, all_z, color="lightgrey", linewidth=1, alpha=0.5, zorder=1)

sc = ax.scatter(xs, ys, zs, c=colors, cmap="RdYlGn_r", s=60, vmin=0, vmax=10, zorder=2)
plt.colorbar(sc, label="ddG (kcal/mol)", shrink=0.6)
ax.set_title("PLP1 Missense Mutations Colored by Structural Destabilization")
ax.set_axis_off()

plt.savefig("../../figures/plp1_mutation_map.png", dpi=300, bbox_inches="tight")
print("Saved figure to figures/plp1_mutation_map.png")