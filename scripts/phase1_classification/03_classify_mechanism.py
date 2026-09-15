import pandas as pd

df = pd.read_csv("../../data/clinvar/plp1_parsed.csv")

pathogenic_labels = ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"]
df_path = df[df["ClinicalSignificance"].isin(pathogenic_labels)].copy()

def assign_mechanism(consequence):
    if consequence == "duplication_cnv":
        return "duplication"
    elif consequence in ["deletion_cnv", "nonsense", "frameshift", "splice_region"]:
        return "loss_of_function"
    elif consequence == "missense":
        return "misfolding_candidate"
    else:
        return "unclassified"

df_path["mechanism"] = df_path["consequence"].apply(assign_mechanism)

df_path.to_csv("../../data/clinvar/plp1_classified.csv", index=False)
print(df_path["mechanism"].value_counts())
print(f"\nTotal pathogenic/likely pathogenic variants: {len(df_path)}")