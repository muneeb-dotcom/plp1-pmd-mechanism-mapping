import subprocess
import glob
import re
import pandas as pd

with open("docking/box_center.txt") as f:
    cx, cy, cz = map(float, f.read().split())

VINA = "tools/vina/vina.exe"
BOX_SIZE = 20

results = []
for lig_path in glob.glob("docking/ligands/*.pdbqt"):
    name = lig_path.split("\\")[-1].replace(".pdbqt", "")
    out_path = f"docking/{name}_out.pdbqt"
    log_path = f"docking/{name}_log.txt"

    cmd = [
        VINA,
        "--receptor", "docking/receptor.pdbqt",
        "--ligand", lig_path,
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(BOX_SIZE), "--size_y", str(BOX_SIZE), "--size_z", str(BOX_SIZE),
        "--out", out_path,
        "--exhaustiveness", "8",
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_path, "w") as f:
        f.write(proc.stdout)

    match = re.search(r"^\s*1\s+(-?\d+\.?\d*)", proc.stdout, re.MULTILINE)
    affinity = float(match.group(1)) if match else None

    results.append({"ligand": name, "affinity_kcal_mol": affinity})
    print(name, "->", affinity, "kcal/mol")

df = pd.DataFrame(results).sort_values("affinity_kcal_mol")
df.to_csv("results/docking_scores.csv", index=False)
print(df.to_string(index=False))