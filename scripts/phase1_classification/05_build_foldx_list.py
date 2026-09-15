import pandas as pd

df = pd.read_csv("../../data/clinvar/plp1_classified.csv")
missense = df[df["mechanism"] == "misfolding_candidate"].dropna(
    subset=["protein_pos", "ref_aa", "alt_aa"]
)

aa3to1 = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E",
    "Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K","Met":"M","Phe":"F",
    "Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V"
}

def build_mutant_string(row):
    ref1 = aa3to1.get(row["ref_aa"])
    alt1 = aa3to1.get(row["alt_aa"])
    if not ref1 or not alt1:
        return None
    # FoldX format: <WT residue><chain><position><mutant residue>;
    return f'{ref1}A{int(row["protein_pos"])}{alt1};'

missense = missense.copy()
missense["foldx_mutation"] = missense.apply(build_mutant_string, axis=1)
clean = missense.dropna(subset=["foldx_mutation"])

clean.to_csv("../../data/structures/missense_variants_clean.csv", index=False)
clean["foldx_mutation"].to_csv("../../data/structures/individual_list.txt", index=False, header=False)

print(f"{len(clean)} missense mutations written to individual_list.txt")