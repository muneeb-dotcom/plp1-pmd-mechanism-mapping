import pandas as pd
import re

df = pd.read_csv("../../data/clinvar/plp1_raw.csv")

def parse_protein_change(name):
    name = str(name)

    match = re.search(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|=)", name)
    if match:
        ref_aa, pos, alt_aa = match.groups()
        pos = int(pos)
        if alt_aa == "Ter":
            consequence = "nonsense"
        elif alt_aa == "=":
            consequence = "synonymous"
        else:
            consequence = "missense"
        return ref_aa, pos, alt_aa, consequence

    fs_match = re.search(r"p\.([A-Za-z]{3})(\d+)fs", name)
    if fs_match:
        ref_aa, pos = fs_match.groups()
        return ref_aa, int(pos), None, "frameshift"

    if re.search(r"x[3-9]\b", name) or "copy number gain" in name.lower():
        return None, None, None, "duplication_cnv"
    if re.search(r"x[01]\b", name) or "copy number loss" in name.lower():
        return None, None, None, "deletion_cnv"
    if re.search(r"c\.\d+[+-]\d+", name):
        return None, None, None, "splice_region"

    return None, None, None, "unknown"

df[["ref_aa", "protein_pos", "alt_aa", "consequence"]] = df["Name"].apply(
    lambda n: pd.Series(parse_protein_change(n))
)

df.to_csv("../../data/clinvar/plp1_parsed.csv", index=False)
print(df["consequence"].value_counts())
print("\nClinicalSignificance sample:")
print(df["ClinicalSignificance"].value_counts(dropna=False).head(10))