import pandas as pd

WT_TOTAL_ENERGY = 1.71  # from RepairPDB's "Total" line on the repaired WT structure

df = pd.read_csv("Average_PLP1_alphafold_Repair.fxout", sep="\t", skiprows=8)
df["mutation_index"] = df["Pdb"].str.extract(r"_(\d+)$").astype(int)
df = df.sort_values("mutation_index").reset_index(drop=True)

mut_df = pd.read_csv("../../data/structures/missense_variants_clean.csv").reset_index(drop=True)
mut_df["mutation_index"] = mut_df.index + 1  # 1-based, matches FoldX's numbering order

merged = mut_df.merge(df[["mutation_index", "total energy"]], on="mutation_index")
merged = merged.rename(columns={"total energy": "mutant_total_energy"})
merged["ddG_kcal_mol"] = merged["mutant_total_energy"] - WT_TOTAL_ENERGY

merged.to_csv("../../data/structures/plp1_structural_disruption_scores.csv", index=False)

print(merged[["Name", "ddG_kcal_mol"]].sort_values("ddG_kcal_mol", ascending=False).head(15))
print(f"\n{len(merged)} mutations scored")
print(f"Destabilizing (ddG > 1.5): {(merged['ddG_kcal_mol'] > 1.5).sum()}")
print(f"Mild (0.5-1.5): {((merged['ddG_kcal_mol'] > 0.5) & (merged['ddG_kcal_mol'] <= 1.5)).sum()}")
print(f"Neutral (<=0.5): {(merged['ddG_kcal_mol'] <= 0.5).sum()}")