import pandas as pd

meta = pd.read_csv("data/heldout/gse96049/nevin_line_metadata.csv")
print(meta.to_string())