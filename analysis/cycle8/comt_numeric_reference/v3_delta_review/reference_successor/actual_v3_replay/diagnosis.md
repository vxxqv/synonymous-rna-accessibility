# First pooled-median mismatch: reference conversion expectation

Diagnosis: the reference expected exact decimal-string rationalization, while the frozen production kernel uses the exact rational value of a Python binary float. The production fraction at this first mismatch is independently reproduced by the latter convention. The reference expectation is incorrect for that frozen executable. This is not a difference in the independently reconstructed contrasts, training membership, site weights, weighted-median tie handling or production fraction serialization.

This is a bounded post-failure diagnosis, not a second observed replay. It recomputed only the one pooled-median statistic for the already-failed code/heavy/fold0 case, from the already-authorized contrasts. It did not rerun normalization, five-fold fitting or scoring, and did not edit either implementation.

## Exact witness

For c.18C>A, assessment fold0 leaves101 training variants at42 sites. Total training weight is42, so half weight is21. Both conversion conventions select exactly the same median boundary: c.321C>G and c.33C>A, with cumulative weight exactly21. The answer is the arithmetic average of those boundary values, with the same tie convention.

| Quantity | Exact value |
|---|---|
| Decimal-string reference | 12281366247441889/20000000000000000 |
| Independently recomputed binary-float convention | 1382758911389497/2251799813685248 |
| Production stored fraction | 1382758911389497/2251799813685248 |
| Absolute exact difference | 7347917049/343597383680000000000000000 |

The exact difference is approximately2.1385253200423905e-17. Both fractions cast to the same binary float,0.6140683123720945. Thus a displayed float comparison would miss this exact rational mismatch; this observation is not a retrospective relaxation of the frozen exact comparison.

The unchanged profiler receipt shows zero discrepancy across all19,337 earlier log/preparation/contrast floating comparisons, including all3,000 compared contrast centers/endpoints across both maps. There is no evidence that the allowed2e-12 raw-arithmetic identity tolerance caused this mismatch: the actual compared floats were identical.

## Executable contract and affected operations

The frozen kernel proxy_models.py (SHA8811af9750e813c798122c7e936b6279a0f3b0e53fea1924cb2d633162a34093), number(x), returns Fraction(x), not Fraction(str(x)). Its validated routine applies this to finite training inputs before weighted medians, PAVA and constrained least squares. Prediction routines convert assessment RNA/load in the same way. The independently authored rational function instead used Fraction(str(x)). The protocol prose says rational arithmetic on supplied stored finite inputs; it does not establish a separate decimal-text rule overriding the pinned executable. The earlier decimal-string interpretation was an erroneous reference assumption.

The five-model production adapter passes the contrast-center Python floats directly into the kernel. serial writes str(Fraction) and float(Fraction), preserving the fitted exact fraction; it does not round the exact result through a float. Evaluation also uses Fraction(observed_T) before subtracting the exact stored prediction. Consequently the discrepancy can affect the expected exact values for pooled/ROI medians, PAVA block means/interpolation, constrained regression statistics/predictions and errors/scores. It is not confined to how a final number is displayed. Site-count weights remain exact reciprocal integers, unchanged by this distinction.

A smallest invented counterexample needs only the value0.1: Fraction(str(0.1))=1/10, whereas Fraction(0.1)=3602879701896397/36028797018963968. Their difference is1/180143985094819840, approximately5.551115123125783e-18, despite equal float displays. A one-value weighted median already distinguishes the two exact-arithmetic expectations. The earlier integer-valued staged examples could not reveal this discrepancy. The general synthetic KKT/partition tests justified the optimization algorithms but did not authenticate a floating-input encoding convention.

## Required next action and limits

Preserve the stopped replay, original seals and this diagnosis. If the independent reference is corrected, make a separately versioned and reviewed post-result correction that converts finite Python floats to their exact binary rational, matching the already-frozen production executable; retain existing exact fractions and raw count-string arithmetic. Add non-dyadic-float synthetic controls for training, prediction and error conversion. Do not change the primary production code, source-derived measurements, models, folds, scientific thresholds or comparison tolerance to obtain agreement.

This diagnosis establishes only the cause of the first stopped comparison. It does not verify all predictions, ranks, losses, coverage or final decision flags. Small error at this one median does not bound every subsequent operation, especially exact-rank branches or ill-conditioned regression. A fresh explicitly authorized complete comparison would still be required after independent review of a documented correction. No computational success, biological gain, mechanism or assay equivalence can be inferred yet.
