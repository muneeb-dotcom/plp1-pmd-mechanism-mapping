import pandas as pd

path = "data/heldout/gse277705/extracted/GSM8527875_Jimpy_perk_Flox_Cre_negative_1_S10.txt.gz"
df = pd.read_csv(path, sep="\t", compression="gzip", nrows=10)
print(df.shape)
print(df.columns.tolist())
print(df.head(10))

full = pd.read_csv(path, sep="\t", compression="gzip")
print(f"\nFull file: {full.shape[0]} rows")