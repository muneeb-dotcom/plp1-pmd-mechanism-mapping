import pandas as pd

de = pd.read_csv("results/plp1_jimpy_deseq_results.csv", index_col=0)
de = de.dropna(subset=["t"])

N = 150
up_genes   = de.nlargest(N, "t").index.tolist()
down_genes = de.nsmallest(N, "t").index.tolist()

with open("results/clue_up_genes.txt", "w") as f:
    f.write("\n".join(up_genes))
with open("results/clue_down_genes.txt", "w") as f:
    f.write("\n".join(down_genes))

print(f"{len(up_genes)} up-regulated genes -> results/clue_up_genes.txt")
print(f"{len(down_genes)} down-regulated genes -> results/clue_down_genes.txt")
print("\nTop 10 up (disease-associated):", up_genes[:10])
print("Top 10 down (disease-associated):", down_genes[:10])

upr_check = ["Eif2ak3", "Hspa5", "Ern1", "Atf4", "Ddit3"]
ranked = de.sort_values("t", ascending=False)
print("\nUPR gene ranks in the disease signature (lower = more disease-upregulated):")
for g in upr_check:
    matches = [i for i, name in enumerate(ranked.index) if str(name).lower() == g.lower()]
    if matches:
        print(f"  {g}: rank {matches[0]+1} of {len(ranked)}")
    else:
        print(f"  {g}: not found in index")