import requests

uniprot_id = "P60201"  # human PLP1
api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

r = requests.get(api_url, timeout=30)
r.raise_for_status()
data = r.json()

pdb_url = data[0]["pdbUrl"]
entry_id = data[0]["entryId"]
print(f"Found model: {entry_id}, downloading from {pdb_url}")

pdb_resp = requests.get(pdb_url, timeout=30)
pdb_resp.raise_for_status()

out_path = "../../data/structures/PLP1_alphafold.pdb"
with open(out_path, "wb") as f:
    f.write(pdb_resp.content)

print(f"Saved to {out_path}, size: {len(pdb_resp.content)} bytes")