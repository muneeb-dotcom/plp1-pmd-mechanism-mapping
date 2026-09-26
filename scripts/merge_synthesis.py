import pandas as pd
import json
import os

lines = []
lines.append("# PLP1/PMD Extended Project — Full Synthesis\n")

# --- Section 1: Part 1 (mechanism map) ---
lines.append("## 1. PLP1 Mechanistic Heterogeneity (Part 1)")
classified = pd.read_csv("data/clinvar/plp1_classified.csv")
strategy = pd.read_csv("results/plp1_mechanism_strategy_map.csv")
lines.append(f"371 pathogenic PLP1 variants classified into 3 mechanism classes: "
             f"{(classified['mechanism']=='loss_of_function').sum()} loss-of-function, "
             f"{(classified['mechanism']=='duplication').sum()} duplication, "
             f"{(classified['mechanism']=='misfolding_candidate').sum()} misfolding candidates. "
             f"FoldX ddG scan confirmed 27 as structurally destabilizing (top: Gly246Trp, ddG=48.3). "
             f"UPR pathway significantly upregulated in a misfolding-mutant model (Wilcoxon p=0.0073). "
             f"FoldX PositionScan identified Leu31Gly as a candidate rescue mutation (-24.3 kcal/mol).\n")

# --- Section 2: Module A (digital twin) ---
lines.append("## 2. Digital Twin (Module A)")
lines.append("Built a mouse oligodendrocyte-lineage reference atlas (Marques et al. 2016, 4,927 QC'd cells), "
             "validated against canonical markers, with a 6-module cell-state framework (myelin_output, "
             "er_stress, opc_identity, apoptosis, lipid_synth) and a maturation pseudotime axis. Modules "
             "validated on independent human data (Jakel et al. 2019). Cross-mechanism expression data "
             "confirmed unavailable for LOF/duplication classes -- twin honestly reports 'no_measured_data' "
             "for those rather than extrapolating. Full repo: "
             "https://github.com/muneeb-dotcom/plp1-pmd-digital-twin\n")

# --- Section 3: Module B (therapeutics) ---
lines.append("## 3. Therapeutic Design (Module B)")
interv = pd.read_csv("module-b-therapeutics/results/intervention_map.csv")
lines.append("Mechanism-matched intervention map built across all 3 classes: AAV replacement (loss-of-"
             "function, dosage-limited MBP promoter to avoid overshoot), ASO knockdown (duplication, "
             "82 candidate 3'UTR target sites), and chaperone/PERK-inhibition (misfolding, target "
             "prioritized via Open Targets, 6 compounds docked against the Gly246Trp/Leu31 pocket, best "
             "-7.7 kcal/mol). CNS-penetrance filtering found no compound fully brain-ready -- a real gap, "
             "not a dead end.\n")

# --- Section 4-7: Module C (already assembled) ---
with open("module-c-validation/results/final_synthesis.md") as f:
    module_c = f.read()
lines.append(module_c)

with open("PROJECT_SYNTHESIS.md", "w") as f:
    f.write("\n".join(lines))

print("Written PROJECT_SYNTHESIS.md")
print(f"Total length: {sum(len(l) for l in lines)} chars")