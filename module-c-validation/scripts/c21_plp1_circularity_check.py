import pandas as pd
import numpy as np
import glob
import os
import re

files = glob.glob("data/heldout/gse96049/extracted/*.fpkm_tracking.gz")

rows = []
for f in files:
    m = re.search(r"GSM\d+_([A-Za-z0-9]+)\.fpkm_tracking", os.path.basename(f))
    sample = m.group(1)
    df = pd.read_csv(f, sep="\t")[["gene_short_name", "FPKM"]]
    plp1_rows = df[df["gene_short_name"] == "PLP1"]
    plp1_fpkm = plp1_rows["FPKM"].mean() if len(plp1_rows) else np.nan
    rows.append({"sample": sample, "PLP1_FPKM": plp1_fpkm})

df = pd.DataFrame(rows).sort_values("sample")
print(df.to_string(index=False))

MYELIN_GENES_NO_PLP1 = ["MBP","MOG","MAG","CNP","MAL","UGT8","MOBP"]
rows2 = []
for f in files:
    m = re.search(r"GSM\d+_([A-Za-z0-9]+)\.fpkm_tracking", os.path.basename(f))
    sample = m.group(1)
    d = pd.read_csv(f, sep="\t")[["gene_short_name", "FPKM"]]
    d.columns = ["gene", "FPKM"]
    d = d.groupby("gene", as_index=False).mean()
    log_fpkm = np.log2(d.set_index("gene")["FPKM"] + 1)
    present = [g for g in MYELIN_GENES_NO_PLP1 if g in log_fpkm.index]
    rows2.append({"sample": sample, "myelin_output_no_plp1": log_fpkm.loc[present].mean()})

df2 = pd.DataFrame(rows2).sort_values("sample")
df2.to_csv("results/plp1_circularity_check.csv", index=False)
print("\nmyelin_output recomputed WITHOUT PLP1 in the panel:")
print(df2.to_string(index=False))