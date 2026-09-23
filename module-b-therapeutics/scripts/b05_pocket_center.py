from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure("mut", "../plp1-pmd-project/tools/foldx/PLP1_alphafold_Repair_12_1.pdb")

pocket_residues = [31, 32, 245, 246, 247]
coords = []
for residue in structure[0]["A"]:
    if residue.id[1] in pocket_residues and "CA" in residue:
        coords.append(residue["CA"].coord)

coords = np.array(coords)
center = coords.mean(axis=0)

print(f"Pocket residues found: {len(coords)} of {len(pocket_residues)}")
print(f"Docking box center: x={center[0]:.2f}, y={center[1]:.2f}, z={center[2]:.2f}")

import os
os.makedirs("docking", exist_ok=True)
with open("docking/box_center.txt", "w") as f:
    f.write(f"{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")