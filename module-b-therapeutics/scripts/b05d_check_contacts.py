from Bio.PDB import PDBParser
import numpy as np
import pandas as pd
import glob

parser = PDBParser(QUIET=True)
structure = parser.get_structure("receptor", "docking/receptor_clean.pdb")

pocket_residues = [31, 32, 245, 246, 247]
pocket_atoms = []
for residue in structure[0]["A"]:
    if residue.id[1] in pocket_residues:
        for atom in residue:
            pocket_atoms.append((residue.id[1], atom.coord))

def parse_top_pose(pdbqt_path):
    coords = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")):
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append(np.array([x, y, z]))
    return coords

rows = []
for out_path in glob.glob("docking/*_out.pdbqt"):
    name = out_path.split("\\")[-1].replace("_out.pdbqt", "")
    lig_coords = parse_top_pose(out_path)

    min_dist = 999
    closest_res = None
    for lig_atom in lig_coords:
        for res_id, pocket_atom in pocket_atoms:
            d = np.linalg.norm(lig_atom - pocket_atom)
            if d < min_dist:
                min_dist = d
                closest_res = res_id

    rows.append({"ligand": name, "min_dist_to_pocket_A": round(min_dist, 2), "closest_residue": closest_res})

df = pd.DataFrame(rows).sort_values("min_dist_to_pocket_A")
scores = pd.read_csv("results/docking_scores.csv")
merged = scores.merge(df, on="ligand")
merged.to_csv("results/docking_scores_with_contacts.csv", index=False)
print(merged.to_string(index=False))