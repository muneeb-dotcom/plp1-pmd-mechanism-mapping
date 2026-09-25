# Module C -- Failure Analysis (Phase C5)

## 1. PLP1-null (PMD12) is the largest outlier despite being clinically mildest
Full PLP1 deletion (PMD12) shows the largest transcriptional displacement from wild-type of all 12 patient lines (2.88, vs range 0.81-2.30 for the rest), even with PLP1 itself excluded from the scoring panel to remove circularity. Clinically this line is the mildest (consistent with the literature: loss-of-function PLP1 mutations are generally milder than toxic gain-of-function point mutations). This means displacement-from-healthy, as scored here, is not a proxy for clinical severity -- it may instead be picking up a distinct axis (e.g. absence of a structural protein vs. presence of a misfolded one triggering ER stress). The two failure modes look transcriptionally similar in this state space but are clinically opposite.

## 2. Severity-rank vs displacement-rank mismatches
Lines ranked by mismatch between clinical severity rank and displacement rank (largest gap first):
       mutation_type  severity  severity_rank  displacement  displacement_rank  rank_gap
10    point_mutation    severe              3      1.086791                2.0      10.0
3      full_deletion      mild              1      2.875384               12.0       8.0
4     point_mutation  moderate              2      0.806033                1.0       7.0
2       triplication    severe              3      1.455526                5.0       7.0
9     point_mutation    severe              3      1.612027                7.0       5.0
0     point_mutation    severe              3      1.782635                8.0       4.0
11  partial_deletion  moderate              2      1.392540                4.0       4.0
1        duplication  moderate              2      1.594189                6.0       2.0
6     point_mutation  moderate              2      1.868523               10.0       2.0
5     point_mutation      mild              1      1.365514                3.0       1.0
7     point_mutation  moderate              2      1.863847                9.0       1.0
8     point_mutation    severe              3      2.296622               11.0       1.0

## 3. Mechanism misclassification pattern (Test 1, binary)
  mechanism_binary   predicted_binary  n
    point_mutation     point_mutation  5
    point_mutation structural_variant  3
structural_variant     point_mutation  4
All 4 structural_variant lines (duplication, triplication, full_deletion, partial_deletion) were misclassified as point_mutation. The state vector has no signal separating these classes -- structural variants are transcriptionally closer to the point-mutation centroid than to each other in this cohort.

## 4. Duplication vs triplication -- dosage non-linearity check
PMD10 (duplication, 2 copies): displacement = 1.594, severity = moderate
PMD11 (triplication, 3 copies): displacement = 1.456, severity = severe
Displacement ratio (3-copy / 2-copy) = 0.91. If dosage effects were linear, tripling copy number relative to duplication (a 1.5x dosage step) would predict a proportionally larger displacement; a ratio far from ~1.5 indicates the model is not capturing dosage response linearly, but n=1 per class makes this a single-pair observation, not a trend.

## Summary
Two of four Part 7 tests failed outright (Test 1, Test 2); the most informative finding is not in the pass/fail table but in Section 1 above: the state vector conflates a mild loss-of-function mechanism with high transcriptional displacement, suggesting displacement-from-healthy and clinical severity are answering different biological questions and should not have been assumed to track each other going in.