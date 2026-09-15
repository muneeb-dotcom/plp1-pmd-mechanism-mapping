import pandas as pd

with open("PS_PLP1_alphafold_Repair_12_1_scanning_output.txt") as f:
    lines = f.readlines()

rows = []
for line in lines:
    parts = line.strip().split("\t")
    if len(parts) == 2:
        mutation, energy = parts
        try:
            rows.append({"mutation": mutation, "total_energy": float(energy)})
        except ValueError:
            continue

df = pd.DataFrame(rows)
df.to_csv("../../results/positionscan_gly246trp_rescue.csv", index=False)

df_sorted = df.sort_values("total_energy")
print("Top 15 most stabilizing rescue candidates near Gly246Trp:\n")
print(df_sorted.head(15).to_string(index=False))

baseline = df[df["mutation"].str.contains("245V")]["total_energy"]
print(f"\nBaseline (unmutated Val245): {baseline.values}")