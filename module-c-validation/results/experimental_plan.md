# Module C -- Experimental Validation Plan (Phase C6)

Three predictions drawn directly from Module C's real findings (not from invented ones), each with a falsification condition specified up front.

## Prediction 1 -- PLP1-null and point-mutant oligodendrocytes fail by different mechanisms
**Basis:** PMD12 (full PLP1 deletion) shows the largest transcriptional displacement of all 12 Nevin lines despite being clinically mildest, even after removing PLP1 itself from the scoring panel (failure_analysis.md, Section 1).
**Hypothesis:** PLP1-null oligodendrocytes fail via loss of a structural/trafficking function, not via the ER-stress/UPR pathway that drives point-mutant pathology.
**Model system:** isogenic iPSC-derived OPCs -- PLP1-null (CRISPR knockout) vs a representative severe point mutant (e.g. PMD1, p.Leu81Arg), same genetic background.
**Readouts:** BiP/CHOP western blot (ER stress), PLP1 subcellular localization by immunofluorescence, MyRF+ maturation counts.
**Sample size:** n=3 independent differentiations per line, matching Nevin et al. 2017 replicate convention.
**Controls:** isogenic wild-type-corrected line for each mutant.
**Falsifies the prediction if:** PLP1-null cells show BiP/CHOP elevation comparable to the point mutant -- would mean both mechanisms converge on ER stress despite the transcriptomic separation seen here.

## Prediction 2 -- PLP1 copy-number response to transcriptional displacement is non-monotonic
**Basis:** triplication (PMD11, severe) shows LOWER displacement (1.456) than duplication (PMD10, moderate, 1.594) -- inverted relative to clinical severity, on n=1 per class (failure_analysis.md, Section 4).
**Hypothesis:** transcriptional displacement does not increase monotonically with PLP1 copy number, unlike clinical severity.
**Model system:** an isogenic CRISPR-engineered PLP1 copy-number allelic series (1, 2, 3, 4 copies) in a single iPSC background, removing the confound of different genetic backgrounds across Nevin's patient lines.
**Readouts:** the same three-module state vector (myelin_output, apoptosis, upr_signaling) used throughout Module C, plus PLP1 protein level by western blot.
**Sample size:** n=3 differentiations per copy-number line (4 lines x 3 = 12 total).
**Controls:** isogenic 1-copy (wild-type) line as baseline.
**Falsifies the prediction if:** displacement increases monotonically with copy number across the full allelic series -- would mean the PMD10/PMD11 inversion was patient-background noise, not a real dosage-response feature.

## Prediction 3 -- PERK inhibition outperforms Ro 25-6981 specifically on myelination, not just survival
**Basis:** PERK inhibition (GSE277705, real data) rescues both survival AND myelin_output together; Ro 25-6981 (Elitt et al. 2018) rescues survival only (myelination_predictions.csv, Test 4). These are different compounds tested in different studies -- never compared head to head.
**Hypothesis:** PERK inhibition rescues the downstream myelination step that Ro 25-6981 fails to rescue, because it acts further downstream in the same ER-stress pathway.
**Model system:** jimpy iPSC-derived OPCs, same in vitro microfiber myelination assay used in Elitt et al. 2018 (Fig 7), head-to-head: vehicle vs Ro 25-6981 vs PERK inhibitor, same plate.
**Readouts:** % microfiber-MBP overlap at day 3 and day 10 (matching Elitt 2018's exact assay), plus DAPI+/MBP+ survival counts as the internal positive control that both compounds should already rescue.
**Sample size:** n=3 replicate wells per compound per timepoint, matching Elitt 2018's own n.
**Controls:** vehicle-treated jimpy, vehicle-treated wild-type.
**Falsifies the prediction if:** PERK inhibitor shows the same survival-without-myelination dissociation as Ro 25-6981 in this side-by-side assay -- would mean the GSE277705 in vivo myelin_output signal doesn't translate to the in vitro fiber assay, and the two datasets aren't actually comparable the way Test 4 assumed.