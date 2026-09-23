import pandas as pd
import json

# --- Step 1: collapse overlapping ASO windows into distinct lead clusters ---
aso = pd.read_csv("../results/aso_targets.csv").sort_values("position").reset_index(drop=True)

clusters = []
current = [aso.iloc[0]]
for i in range(1, len(aso)):
    if aso.iloc[i]["position"] - current[-1]["position"] <= 3:
        current.append(aso.iloc[i])
    else:
        clusters.append(current)
        current = [aso.iloc[i]]
clusters.append(current)

lead_sites = []
for cluster in clusters:
    best = max(cluster, key=lambda r: r["accessibility"])
    lead_sites.append({
        "position": best["position"],
        "sequence_dna": best["sequence_dna"],
        "gc_content": best["gc_content"],
        "accessibility": best["accessibility"],
        "cluster_size": len(cluster),
    })

lead_df = pd.DataFrame(lead_sites).sort_values("accessibility", ascending=False)
lead_df.to_csv("../results/aso_lead_sites.csv", index=False)

print(f"{len(clusters)} distinct ASO target regions identified from {len(aso)} passing windows")
print("\nTop 5 lead sites (one per genomic cluster):")
print(lead_df.head(5).to_string(index=False))

# --- Step 2: build the final mechanism-matched intervention map ---
with open("data/correction_targets.json") as f:
    targets = json.load(f)

target_priority = pd.read_csv("../results/target_priority.csv")
cns = pd.read_csv("../results/cns_filtered_candidates.csv")

rows = []
rows.append({
    "mechanism": "loss_of_function", "n_variants": targets["loss_of_function"]["n_variants"],
    "modality": "AAV gene replacement (dosage-limited promoter, see Part 4 design notes)",
    "lead_candidates": "n/a — cassette design, not a compound list",
    "evidence_tier": "established clinical modality (CNS AAV gene therapy precedent)",
})
rows.append({
    "mechanism": "duplication", "n_variants": targets["duplication"]["n_variants"],
    "modality": "ASO, partial knockdown, 3'UTR target (shared by PLP1/DM20)",
    "lead_candidates": f"{len(lead_df)} candidate sites; top: pos {lead_df.iloc[0]['position']} ({lead_df.iloc[0]['sequence_dna']})",
    "evidence_tier": "strong — PLP1 ASO suppression already validated (Elitt et al. 2020, Nature)",
})
best_tractable = target_priority[target_priority["tractable_smallmolecule"]]["symbol"].tolist()
rows.append({
    "mechanism": "misfolding_confirmed", "n_variants": targets["misfolding_confirmed"]["n_variants"],
    "modality": f"Chaperone/UPR modulator targeting {', '.join(best_tractable)}",
    "lead_candidates": f"{len(cns)} compounds scored; best tier: {cns['cns_tier'].mode()[0]} (none fully CNS-favourable — needs med-chem optimization)",
    "evidence_tier": "preclinical (ISRIB-class ISR modulators in development)",
})
rows.append({
    "mechanism": "misfolding_uncertain", "n_variants": targets["misfolding_uncertain"]["n_variants"],
    "modality": "undetermined",
    "lead_candidates": "none — requires functional validation first",
    "evidence_tier": "n/a",
})

intervention_map = pd.DataFrame(rows)
intervention_map.to_csv("../results/intervention_map.csv", index=False)

print("\n" + "=" * 70)
print("FINAL MECHANISM-MATCHED INTERVENTION MAP")
print("=" * 70)
for _, r in intervention_map.iterrows():
    print(f"\n{r['mechanism']}  (n={r['n_variants']})")
    print(f"  modality:   {r['modality']}")
    print(f"  candidates: {r['lead_candidates']}")
    print(f"  evidence:   {r['evidence_tier']}")