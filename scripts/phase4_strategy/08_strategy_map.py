import pandas as pd

classified = pd.read_csv("../../data/clinvar/plp1_classified.csv")
structural = pd.read_csv("../../data/structures/plp1_structural_disruption_scores.csv")

structural_lookup = structural.set_index("Name")["ddG_kcal_mol"].to_dict()
classified["ddG_kcal_mol"] = classified["Name"].map(structural_lookup)

def refine_mechanism(row):
    if row["mechanism"] != "misfolding_candidate":
        return row["mechanism"]
    ddg = row["ddG_kcal_mol"]
    if pd.isna(ddg):
        return "misfolding_candidate_unscored"
    elif ddg > 1.5:
        return "misfolding_confirmed"
    else:
        return "misfolding_uncertain"

classified["mechanism_refined"] = classified.apply(refine_mechanism, axis=1)

strategy_map = {
    "duplication": "Suppression therapy (ASO/RNAi knockdown of PLP1 mRNA)",
    "loss_of_function": "Replacement therapy (AAV-mediated PLP1 gene delivery)",
    "misfolding_confirmed": "Stabilization therapy (chaperone/small-molecule folding correctors or UPR-modulators, e.g. ISRIB-class compounds)",
    "misfolding_uncertain": "Uncertain — pathogenic but structurally mild; needs functional validation before strategy assignment",
    "misfolding_candidate_unscored": "Not structurally scored (excluded from Phase 2 FoldX batch)",
    "unclassified": "Insufficient annotation to assign mechanism"
}

classified["therapeutic_strategy"] = classified["mechanism_refined"].map(strategy_map)

classified.to_csv("../../results/plp1_mechanism_strategy_map.csv", index=False)

print("Refined mechanism counts:")
print(classified["mechanism_refined"].value_counts())
print("\nTherapeutic strategy breakdown:")
print(classified["therapeutic_strategy"].value_counts())