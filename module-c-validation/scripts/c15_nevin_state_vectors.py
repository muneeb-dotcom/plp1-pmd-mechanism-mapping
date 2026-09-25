import pandas as pd
import numpy as np
import glob
import os
import re

MODULES = {
    "myelin_output": ["MBP","MOG","MAG","CNP","MAL","PLP1","UGT8","MOBP"],
    "apoptosis": ["CASP3","BAX","BCL2L11","TP53","CDKN1A"],
    "upr_signaling": ["ATF4","DDIT3","EIF2AK3","XBP1","ATF6"],
}

files = glob.glob("data/heldout/gse96049/extracted/*.fpkm_tracking.gz")
print(f"{len(files)} files found")

rows = []
for f in files:
    m = re.search(r"GSM\d+_([A-Za-z0-9]+)\.fpkm_tracking", os.path.basename(f))
    sample = m.group(1)
    df = pd.read_csv(f, sep="\t")[["gene_short_name", "FPKM"]]
    df.columns = ["gene", "FPKM"]
    df = df.groupby("gene", as_index=False).mean()
    log_fpkm = np.log2(df.set_index("gene")["FPKM"] + 1)

    row = {"sample": sample}
    for module, genes in MODULES.items():
        present = [g for g in genes if g in log_fpkm.index]
        row[module] = log_fpkm.loc[present].mean() if present else np.nan
        row[f"{module}_ngenes"] = len(present)
    rows.append(row)

state = pd.DataFrame(rows).sort_values("sample")
state.to_csv("results/nevin_state_vectors.csv", index=False)
print(state.to_string(index=False))