\# PLP1/PMD Twin Validation (Module C)



Tests the Module A digital twin on real held-out data it was never

fitted on.



Part of the extended PLP1/PMD project:

\- Part 1: https://github.com/muneeb-dotcom/plp1-pmd-mechanism-mapping

\- Module A: https://github.com/muneeb-dotcom/plp1-pmd-digital-twin

\- Module B: (therapeutic design)



\## Key finding

Original planned validation sets (Nevin et al. 2017, Elitt et al.

2020\) have no public transcriptomic data. Found and used GSE277705

instead (Chen et al. 2025) - real jimpy mice with PERK genetically

knocked out, directly testing Module B's top-priority target.



Result: myelin\_output increased strongly (+1.0 log2CPM) under PERK

knockout, matching the source paper. Original er\_stress module was

found to conflate UPR signaling (correctly decreased) with chaperone

capacity (compensatorily increased) - module was split and re-validated.

See results/split\_module\_validation.csv and results/perk\_knockout\_validation.csv.



\## Structure

\- `scripts/` — c00 through c07

\- `results/` — validation outputs

\- `data/gene\_symbol\_to\_ensembl.json` — verified mouse gene ID mapping



\## Reproduce

Raw GSE277705 files gitignored (\~12MB, regenerable via scripts/c01).

