import pandas as pd

df = pd.read_csv("../../data/structures/plp1_structural_disruption_scores.csv")
target = df.sort_values("ddG_kcal_mol", ascending=False).iloc[0]

print(target[["Name", "protein_pos", "ref_aa", "alt_aa", "ddG_kcal_mol"]])