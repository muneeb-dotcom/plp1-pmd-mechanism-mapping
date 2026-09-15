from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure("mut", "PLP1_alphafold_Repair_12_1.pdb")

target_ca = None
for residue in structure[0]["A"]:
    if residue.id[1] == 246:
        target_ca = residue["CA"].coord
        break

aa3to1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

neighbors = []
for residue in structure[0]["A"]:
    if "CA" not in residue or residue.id[1] == 246:
        continue
    dist = np.linalg.norm(residue["CA"].coord - target_ca)
    if dist <= 8.0:
        aa1 = aa3to1.get(residue.get_resname())
        if aa1:
            neighbors.append((residue.id[1], aa1, dist))

neighbors.sort(key=lambda x: x[2])
print(f"{len(neighbors)} residues within 8A of position 246:\n")
for pos, aa1, dist in neighbors:
    print(f"{aa1}A{pos}  ({dist:.1f} A away)")

position_string = ",".join(f"{aa1}A{pos}a" for pos, aa1, dist in neighbors)
print("\nFoldX PositionScan string:")
print(position_string)