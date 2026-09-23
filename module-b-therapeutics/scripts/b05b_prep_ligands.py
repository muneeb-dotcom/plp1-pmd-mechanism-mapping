import os
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess

CANDIDATES = [
    {"name": "PERK_CHEMBL2171126", "smiles": "Cn1cc(-c2ccc3c(c2)CCN3C(=O)Cc2cc(F)cc(F)c2F)c2c(N)ncnc21"},
    {"name": "PERK_CHEMBL2171125", "smiles": "Cn1cc(-c2ccc3c(c2)CCN3C(=O)Cc2cc(F)cc(C(F)(F)F)c2)c2c(N)ncnc21"},
    {"name": "PERK_CHEMBL1667910", "smiles": "Oc1cccc2cc(Nc3c(-c4ncccn4)oc4cnccc34)ccc12"},
    {"name": "Lestaurtinib", "smiles": "C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)CNC4=O"},
    {"name": "Staurosporine", "smiles": "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4"},
    {"name": "IRE1_CHEMBL3356002", "smiles": "Cc1cc(Nc2nc3ccccc3[nH]2)c2ccccc2c1Oc1ncccc1-c1ccnc(N[C@H]2CCCNC2)n1"},
]

os.makedirs("docking/ligands", exist_ok=True)

for c in CANDIDATES:
    mol = Chem.MolFromSmiles(c["smiles"])
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    pdb_path = f"docking/ligands/{c['name']}.pdb"
    pdbqt_path = f"docking/ligands/{c['name']}.pdbqt"
    Chem.MolToPDBFile(mol, pdb_path)

    subprocess.run([
        "obabel", pdb_path, "-O", pdbqt_path,
        "--partialcharge", "gasteiger"
    ], check=True)

    print("Prepared", pdbqt_path)