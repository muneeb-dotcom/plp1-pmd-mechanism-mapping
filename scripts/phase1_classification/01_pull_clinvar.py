from Bio import Entrez
import pandas as pd
import json
import time

Entrez.email = "shahzadmuneeb15@gmail.com"

search = Entrez.esearch(db="clinvar", term="PLP1[gene]", retmax=5000)
record = Entrez.read(search)
ids = record["IdList"]
print(f"Found {len(ids)} ClinVar records for PLP1")

rows = []
batch_size = 50
for i in range(0, len(ids), batch_size):
    batch = ids[i:i + batch_size]
    handle = Entrez.esummary(db="clinvar", id=",".join(batch), retmode="json")
    raw_text = handle.read()
    try:
        data = json.loads(raw_text)
    except json.JSONDecodeError:
        print("Batch", i, "returned non-JSON:", raw_text[:300])
        time.sleep(2)
        continue

    if "result" not in data:
        print("Batch", i, "error response:", data)
        time.sleep(2)
        continue

    for uid in data["result"]["uids"]:
        rec = data["result"][uid]
        cs = rec.get("germline_classification") or rec.get("clinical_significance") or {}
        traits = rec.get("trait_set", []) or []
        phenotypes = ";".join(t.get("trait_name", "") for t in traits)

        chrom, start, stop = None, None, None
        for vs in rec.get("variation_set", []) or []:
            for loc in vs.get("variation_loc", []) or []:
                if loc.get("assembly_name") == "GRCh38":
                    chrom = loc.get("chr")
                    start = loc.get("start")
                    stop = loc.get("stop")

        rows.append({
            "Name": rec.get("title"),
            "Type": rec.get("obj_type"),
            "ClinicalSignificance": cs.get("description"),
            "PhenotypeList": phenotypes,
            "ChromosomeAccession": chrom,
            "Start": start,
            "Stop": stop,
        })
    print(f"Processed batch {i}-{i+len(batch)}, total rows so far: {len(rows)}")
    time.sleep(0.5)

df = pd.DataFrame(rows).drop_duplicates()
df.to_csv("../../data/clinvar/plp1_raw.csv", index=False)
print(len(df), "PLP1 variants pulled")