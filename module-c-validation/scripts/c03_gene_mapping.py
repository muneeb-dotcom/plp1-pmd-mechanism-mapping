import requests
import json
import time

with open("../plp1-twin/results/state_modules.json") as f:
    ref = json.load(f)

all_genes = set()
for module, genes in ref["genes"].items():
    all_genes.update(genes)
all_genes = sorted(all_genes)
print(f"{len(all_genes)} unique module genes to map")

mapping = {}
failed = []

for gene in all_genes:
    url = f"https://rest.ensembl.org/xrefs/symbol/mus_musculus/{gene}?content-type=application/json"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        hits = r.json()
    except Exception as e:
        failed.append((gene, str(e)))
        continue

    gene_hits = [h for h in hits if h.get("id", "").startswith("ENSMUSG")]
    if gene_hits:
        mapping[gene] = gene_hits[0]["id"]
    else:
        failed.append((gene, "no ENSMUSG hit"))
    time.sleep(0.1)

print(f"\nMapped {len(mapping)} of {len(all_genes)} genes")
if failed:
    print("Failed:", failed)

with open("data/gene_symbol_to_ensembl.json", "w") as f:
    json.dump(mapping, f, indent=2)