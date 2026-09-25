import pandas as pd
import numpy as np
import json

with open("../plp1-twin/results/state_modules.json") as f:
    ref = json.load(f)
MODULE_GENES = ref["genes"]

files = {
    "WT10": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034716_WT10-OPC-2.fpkm_tracking.gz",
    "WT13": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034717_WT13-OPC-2.fpkm_tracking.gz",
    "WT14": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034718_WT14-OPC-3.fpkm_tracking.gz",
    "jp14": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034719_jp14-OPC-3.fpkm_tracking.gz",
    "jp15": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034720_jp15-OPC-1.fpkm_tracking.gz",
    "jp16": "../plp1-pmd-project/scripts/GSE111605/extracted/GSM3034721_jp16-OPC-2.fpkm_tracking.gz",
}

def read_fpkm(path, name):
    df = pd.read_csv(path, sep="\t")[["gene_short_name", "FPKM"]]
    df.columns = ["gene", name]
    return df.groupby("gene", as_index=False).mean()

dfs = [read_fpkm(p, n) for n, p in files.items()]
expr = dfs[0]
for d in dfs[1:]:
    expr = expr.merge(d, on="gene", how="outer")
expr = expr.set_index("gene").fillna(0)
log_expr = np.log2(expr + 1)

wt_cols = ["WT10", "WT13", "WT14"]
jp_cols = ["jp14", "jp15", "jp16"]

print("HELD-OUT TIMEPOINT: OPC stage (D1-T3 timepoint was used to build the twin; this was never touched)\n")
rows = []
for module, genes in MODULE_GENES.items():
    present = [g for g in genes if g in log_expr.index]
    if not present:
        continue
    wt_mean = log_expr.loc[present, wt_cols].mean(axis=0).mean()
    jp_mean = log_expr.loc[present, jp_cols].mean(axis=0).mean()
    rows.append({"module": module, "wt_mean": wt_mean, "jimpy_mean": jp_mean, "delta": jp_mean - wt_mean})

df = pd.DataFrame(rows)
print(df.round(3).to_string(index=False))

d1t3 = pd.read_csv("../plp1-twin/results/pmd_state_vectors.csv")
print("\nCompare direction against original D1-T3 result (used to build the twin):")
merged = df.merge(d1t3[["module", "delta_log2fpkm"]], on="module")
merged["same_direction"] = np.sign(merged["delta"]) == np.sign(merged["delta_log2fpkm"])
merged.columns = ["module", "opc_wt", "opc_jimpy", "opc_delta", "d1t3_delta", "same_direction"]
print(merged.round(3).to_string(index=False))
print(f"\n{merged['same_direction'].sum()} of {len(merged)} modules replicate direction at the held-out OPC timepoint")

merged.to_csv("results/opc_timepoint_replication.csv", index=False)