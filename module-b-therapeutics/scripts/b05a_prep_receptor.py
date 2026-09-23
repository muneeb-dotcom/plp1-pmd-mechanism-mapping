from Bio.PDB import PDBParser, PDBIO, Select
import subprocess

receptor_pdb = "../plp1-pmd-project/tools/foldx/PLP1_alphafold_Repair_12_1.pdb"
clean_pdb = "docking/receptor_clean.pdb"
out_pdbqt = "docking/receptor.pdbqt"

class ProteinOnly(Select):
    def accept_residue(self, residue):
        return residue.id[0] == " "  # keep only standard amino acid residues

parser = PDBParser(QUIET=True)
structure = parser.get_structure("receptor", receptor_pdb)

io = PDBIO()
io.set_structure(structure)
io.save(clean_pdb, ProteinOnly())

print("Clean PDB written:", clean_pdb)

subprocess.run([
    "obabel", clean_pdb, "-O", out_pdbqt,
    "-xr", "--partialcharge", "gasteiger", "-h"
], check=True)

print("Receptor written to", out_pdbqt)