import urllib.request
import os

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96049/suppl/GSE96049_RAW.tar"
outdir = "data/heldout/gse96049"
os.makedirs(outdir, exist_ok=True)
out_path = os.path.join(outdir, "GSE96049_RAW.tar")

print("Downloading", url)
urllib.request.urlretrieve(url, out_path)
print("Saved to", out_path, "-", os.path.getsize(out_path), "bytes")