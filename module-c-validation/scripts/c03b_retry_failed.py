import requests
import json
import time

with open("data/gene_symbol_to_ensembl.json") as f:
    mapping = json.load(f)

failed_genes = ["Cyp51", "Derl1", "Dnajb9", "Edem1", "Fdft1", "Herpud1", "Hmgcr",
                "Hspa5", "Hyou1", "Idi1", "Mag", "Mal", "Manf", "Plp1", "Sel1l",
                "Sox10", "Sqle", "Trp53", "Xbp1"]

still_failed = []
for gene in failed_genes:
    url = f"https://rest.ensembl.org/xrefs/symbol/mus_musculus/{gene}?content-type=application/json"
    success = False
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                hits = r.json()
                gene_hits = [h for h in hits if h.get("id", "").startswith("ENSMUSG")]
                if gene_hits:
                    mapping[gene] = gene_hits[0]["id"]
                    success = True
                break
        except Exception:
            pass
        time.sleep(2)
    if not success:
        still_failed.append(gene)
    time.sleep(0.5)

print(f"Recovered {len(failed_genes) - len(still_failed)} of {len(failed_genes)}")
print("Still failed:", still_failed)
print(f"\nTotal mapped now: {len(mapping)} of 40")

with open("data/gene_symbol_to_ensembl.json", "w") as f:
    json.dump(mapping, f, indent=2)