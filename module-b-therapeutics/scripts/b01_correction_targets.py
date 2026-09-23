import json

# Correction targets per mechanism class.
# "direction" and "band" define what "corrected" means numerically —
# this stops later phases from silently assuming "more is better".
CORRECTION_TARGETS = {
    "loss_of_function": {
        "n_variants": 173,
        "disease_state": "PLP1 absent or non-functional; axonal degeneration without demyelination",
        "target_state": "Restore PLP1 expression toward physiological level",
        "direction": "increase",
        "dosage_band_pct_of_normal": [80, 120],
        "overshoot_risk": "Exceeding the band recreates the duplication phenotype (Karim et al. 2007, Glia)",
    },
    "duplication": {
        "n_variants": 63,
        "disease_state": "PLP1 overexpressed; protein accumulates in endolysosomes/autophagic vacuoles; MBP reduced",
        "target_state": "Reduce PLP1 dosage toward physiological level",
        "direction": "decrease",
        "dosage_band_pct_of_normal": [80, 120],
        "overshoot_risk": "Over-suppression recreates the loss-of-function phenotype",
    },
    "misfolding_confirmed": {
        "n_variants": 27,
        "disease_state": "PLP1 expressed but ER-retained; UPR/ER-stress active; oligodendrocyte death",
        "target_state": "Improve folding/trafficking efficiency or dampen chronic UPR signalling",
        "direction": "redirect",
        "dosage_band_pct_of_normal": None,
        "overshoot_risk": "Systemic UPR inhibition is dangerous — the UPR is essential in secretory tissues",
    },
    "misfolding_uncertain": {
        "n_variants": 56,
        "disease_state": "Pathogenic but structurally mild (ddG <= 1.5) — mechanism not resolved by folding-stability alone",
        "target_state": "Undetermined — needs functional validation before a target can be defined",
        "direction": None,
        "dosage_band_pct_of_normal": None,
        "overshoot_risk": None,
    },
}

with open("data/correction_targets.json", "w") as f:
    json.dump(CORRECTION_TARGETS, f, indent=2)

print("Correction targets written for", len(CORRECTION_TARGETS), "mechanism classes")
for mech, info in CORRECTION_TARGETS.items():
    print(f"\n{mech} (n={info['n_variants']})")
    print(f"  direction: {info['direction']}")
    print(f"  band: {info['dosage_band_pct_of_normal']}")