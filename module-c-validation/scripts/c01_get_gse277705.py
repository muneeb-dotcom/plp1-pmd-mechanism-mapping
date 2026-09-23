import requests
import tarfile
import os

os.makedirs("data/heldout/gse277705", exist_ok=True)
tar_path = "data/heldout/gse277705/GSE277705_RAW.tar"
url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277705/suppl/GSE277705_RAW.tar"
headers = {"User-Agent": "Mozilla/5.0"}

r = requests.head(url, headers=headers, timeout=30)
total_size = int(r.headers.get("content-length", 0))

for attempt in range(1, 11):
    existing = os.path.getsize(tar_path) if os.path.exists(tar_path) else 0
    if existing >= total_size:
        break
    print(f"Attempt {attempt}: resuming from {existing/(1024*1024):.1f} MB / {total_size/(1024*1024):.1f} MB")
    range_headers = {**headers, "Range": f"bytes={existing}-"}
    try:
        with requests.get(url, headers=range_headers, stream=True, timeout=60) as resp:
            mode = "ab" if existing > 0 else "wb"
            with open(tar_path, mode) as f:
                for chunk in resp.iter_content(256 * 1024):
                    f.write(chunk)
    except Exception as e:
        print(f"Dropped: {e}, retrying...")
        continue

final = os.path.getsize(tar_path)
print(f"Final size: {final/(1024*1024):.1f} MB / {total_size/(1024*1024):.1f} MB")

if final >= total_size:
    with tarfile.open(tar_path) as tar:
        members = tar.getnames()
        print(f"{len(members)} files in archive")
        tar.extractall("data/heldout/gse277705/extracted")
    jimpy_files = [m for m in members if "Jimpy_perk" in m]
    print(f"\n{len(jimpy_files)} Jimpy_perk arm files:")
    for f in jimpy_files:
        print(" ", f)
else:
    print("INCOMPLETE -- run this script again to keep resuming.")