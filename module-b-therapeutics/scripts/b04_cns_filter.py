import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

# Real compounds pulled from ChEMBL bioactivity data (Phase B3/B4).
# Each entry cites its ChEMBL molecule ID and measured potency so the
# source is traceable, not just a SMILES string from nowhere.
CANDIDATES = [
    # --- PERK (EIF2AK3) inhibitors, target CHEMBL6030 ---
    {"name": "GSK-family PERK inhibitor (CHEMBL2171126)", "target": "EIF2AK3",
     "chembl_id": "CHEMBL2171126", "potency_nM": 0.2, "assay": "IC50",
     "smiles": "Cn1cc(-c2ccc3c(c2)CCN3C(=O)Cc2cc(F)cc(F)c2F)c2c(N)ncnc21"},
    {"name": "PERK inhibitor analog (CHEMBL2171125)", "target": "EIF2AK3",
     "chembl_id": "CHEMBL2171125", "potency_nM": 0.3, "assay": "IC50",
     "smiles": "Cn1cc(-c2ccc3c(c2)CCN3C(=O)Cc2cc(F)cc(C(F)(F)F)c2)c2c(N)ncnc21"},
    {"name": "Early PERK lead (CHEMBL1667910)", "target": "EIF2AK3",
     "chembl_id": "CHEMBL1667910", "potency_nM": 4.1, "assay": "IC50",
     "smiles": "Oc1cccc2cc(Nc3c(-c4ncccn4)oc4cnccc34)ccc12"},
    # --- IRE1 (ERN1) inhibitors, target CHEMBL1163101 ---
    {"name": "Lestaurtinib", "target": "ERN1",
     "chembl_id": "CHEMBL603469", "potency_nM": 5.7, "assay": "Kd",
     "smiles": "C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)CNC4=O"},
    {"name": "Staurosporine", "target": "ERN1",
     "chembl_id": "CHEMBL388978", "potency_nM": 21.0, "assay": "Kd",
     "smiles": "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4"},
    {"name": "Selective IRE1 inhibitor (CHEMBL3356002)", "target": "ERN1",
     "chembl_id": "CHEMBL3356002", "potency_nM": 5.1, "assay": "IC50",
     "smiles": "Cc1cc(Nc2nc3ccccc3[nH]2)c2ccccc2c1Oc1ncccc1-c1ccnc(N[C@H]2CCCNC2)n1"},
]

def cns_mpo_descriptors(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return {
        "MW":   round(Descriptors.MolWt(m), 1),
        "logP": round(Crippen.MolLogP(m), 2),
        "TPSA": round(Descriptors.TPSA(m), 1),
        "HBD":  Descriptors.NumHDonors(m),
        "rotB": Descriptors.NumRotatableBonds(m),
    }

def cns_tier(d):
    """Wager et al. CNS MPO-style desirability flags.
    Favourable ranges: MW<=360, logP<=3, TPSA 40-90, HBD<=1, rotB<=5."""
    flags = {
        "MW_ok":   d["MW"] <= 360,
        "logP_ok": d["logP"] <= 3,
        "TPSA_ok": 40 <= d["TPSA"] <= 90,
        "HBD_ok":  d["HBD"] <= 1,
        "rotB_ok": d["rotB"] <= 5,
    }
    n_pass = sum(flags.values())
    if n_pass >= 4:
        return "CNS-favourable", flags
    elif n_pass >= 2:
        return "borderline", flags
    else:
        return "CNS-unfavourable", flags

rows = []
for c in CANDIDATES:
    desc = cns_mpo_descriptors(c["smiles"])
    if desc is None:
        continue
    tier, flags = cns_tier(desc)
    failed = [k.replace("_ok", "") for k, v in flags.items() if not v]
    rows.append({
        **c, **desc,
        "cns_tier": tier,
        "failed_criteria": ", ".join(failed) if failed else "none",
    })

df = pd.DataFrame(rows)
df.to_csv("../results/cns_filtered_candidates.csv", index=False)

print(df[["name", "target", "potency_nM", "MW", "logP", "TPSA", "cns_tier", "failed_criteria"]]
      .to_string(index=False))

print("\nTier counts:")
print(df["cns_tier"].value_counts())