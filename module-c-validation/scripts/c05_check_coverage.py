import json

with open("../plp1-twin/results/state_modules.json") as f:
    ref = json.load(f)
with open("data/gene_symbol_to_ensembl.json") as f:
    sym2ens = json.load(f)

import pandas as pd
mat_genes = set(pd.read_csv("data/heldout/gse277705/extracted/GSM8527875_Jimpy_perk_Flox_Cre_negative_1_S10.txt.gz",
                            sep="\t", compression="gzip", index_col="Gene").index)

for module, genes in ref["genes"].items():
    print(f"\n{module}:")
    for g in genes:
        ens = sym2ens.get(g, "NOT_MAPPED")
        present = ens in mat_genes if ens != "NOT_MAPPED" else False
        print(f"  {g} -> {ens} -> {'FOUND' if present else 'missing'}")