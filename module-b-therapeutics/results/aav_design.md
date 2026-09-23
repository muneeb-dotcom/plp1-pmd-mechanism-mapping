\# AAV Replacement Design — Loss-of-Function PLP1 (n=173)



\## Cassette size

| Element | Size |

|---|---|

| PLP1 CDS (NM\_000533.5, confirmed via GenBank) | 834 bp |

| Promoter (cell-type-specific, see below) | \~1,200 bp |

| polyA signal | \~250 bp |

| ITRs (x2) | \~300 bp |

| \*\*Total\*\* | \*\*\~2,584 bp\*\* — well within AAV's \~4.7 kb limit |



Packaging is not a constraint for this gene.



\## Capsid

| Capsid | Tropism | Route | Clinical precedent |

|---|---|---|---|

| AAV9 | Crosses BBB, broad CNS transduction | IV or intrathecal | Used in Zolgensma (SMA) |

| AAV.PHP.eB | Very strong CNS transduction | IV | Characterized mainly in specific mouse strains — translation to human uncertain |

| AAVrh10 | Good CNS/oligodendrocyte transduction | Intrathecal/ICV | Used in several CNS gene therapy trials |



\*\*Recommendation: AAV9\*\*, on existing clinical precedent, over AAV.PHP.eB's stronger but less-validated-in-human tropism.



\## Promoter — the critical choice

Karim et al. 2007 (Glia) showed \*Plp1\*-transgenic mice: hemizygous (mild excess) tolerate the increase; homozygous (larger excess) develop hypomyelination and reduced MBP — i.e., \*\*overshoot recreates the duplication phenotype\*\*.



| Promoter | Strength | Risk |

|---|---|---|

| CAG/CMV | Very strong, ubiquitous | High overshoot risk |

| MBP/CNP | Oligodendrocyte-specific, moderate | Lower overshoot risk |

| Endogenous PLP1 promoter | Native regulation | Likely lowest overshoot risk, but weaker/less characterized in AAV context |



\*\*Recommendation: MBP promoter.\*\* Cell-type-restricted and moderate-strength, directly targeting the dosage band defined in Phase B1 (80–120% of normal), rather than a strong ubiquitous promoter that risks recreating the duplication phenotype (Table in Phase B1: `data/correction\_targets.json`).



\## Predicted molecular response

Dose–response is a \*\*band, not a point target\*\* (per B1). Restoring PLP1 to 80–120% of physiological expression is the target window; both under- and over-correction have documented phenotypic consequences (LOF and duplication respectively).



\## Stated limitations

This is a proposed design rationale, not a validated construct. No expression-level, immunogenicity, biodistribution, or tropism data specific to this cassette exist. AAV.PHP.eB's translation risk and the MBP promoter's AAV-context strength are both explicitly flagged as open questions requiring empirical testing.

