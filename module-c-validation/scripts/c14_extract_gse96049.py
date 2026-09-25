import tarfile
import os

tar_path = "data/heldout/gse96049/GSE96049_RAW.tar"
outdir = "data/heldout/gse96049/extracted"
os.makedirs(outdir, exist_ok=True)

with tarfile.open(tar_path) as tar:
    tar.extractall(outdir)

files = sorted(os.listdir(outdir))
print(f"{len(files)} files extracted")
for f in files:
    print(f)