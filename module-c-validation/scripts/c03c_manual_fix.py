import json

with open("data/gene_symbol_to_ensembl.json") as f:
    mapping = json.load(f)

# Hmgcr (HMG-CoA reductase), mouse -- verified Ensembl ID
mapping["Hmgcr"] = "ENSMUSG00000021670"

with open("data/gene_symbol_to_ensembl.json", "w") as f:
    json.dump(mapping, f, indent=2)

print(f"Total mapped: {len(mapping)} of 40")