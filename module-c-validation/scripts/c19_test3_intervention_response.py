import pandas as pd

rows = [{
    "test": "intervention_response",
    "reference": "Elitt et al. 2020, Nature, PMID 32610343",
    "intervention": "Plp1-targeting ASO, single ICV dose, postnatal jimpy mice",
    "reported_direction": "oligodendrocyte numbers up, myelination up, motor/respiratory improved, lifespan extended -- survival AND myelination co-improve",
    "predicted_by_project": "not computable -- no fitted twin exists in this pipeline that outputs a scalar prediction for PLP1-ASO; no RNA-seq deposited for this arm to score against",
    "metric_spec": "directional concordance (per module C4 Test 3 spec)",
    "status": "not_measurable",
}]

df = pd.DataFrame(rows)
df.to_csv("results/test3_intervention_response.csv", index=False)
print(df.to_string(index=False))