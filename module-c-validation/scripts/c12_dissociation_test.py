import pandas as pd

published = {
    "source": "Elitt et al. 2018, Stem Cell Reports, PMID 30146490 (Fig 1B, Fig 7D-E)",
    "note": "Ro 25-6981 arm is immunocytochemistry, not RNA-seq -- no file exists to download",
    "survival_rescued": True,
    "myelination_rescued": False,
}

perk = pd.read_csv("results/gate_math_validation.csv", index_col=0)
survival_up = perk.loc["KO", "apoptosis"] < perk.loc["active", "apoptosis"]
myelin_up = perk.loc["KO", "myelin_output"] > perk.loc["active", "myelin_output"]

rows = [{
    "test": "dissociation_reproduction",
    "reference": published["source"],
    "reference_note": published["note"],
    "reference_pattern": "survival up, myelination NOT up",
    "tested_intervention": "PERK inhibition (GSE277705) -- different compound than Ro 25-6981",
    "survival_up_in_data": survival_up,
    "myelination_up_in_data": myelin_up,
    "matches_published_dissociation": bool(survival_up and not myelin_up),
    "result": "PERK inhibition rescues BOTH survival and myelination -- does not reproduce the Ro 25-6981 dissociation. Different mechanism, different real outcome.",
}]

df = pd.DataFrame(rows)
df.to_csv("results/myelination_predictions.csv", index=False)
print(df.to_string(index=False))