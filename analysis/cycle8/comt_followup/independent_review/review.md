# Independent exploratory COMT diagnostic review

Verdict: PASS for the stated finite fixed-prediction computations. No substantive numerical mismatch, identity loss, normalization change or dropped target was found. This does not establish biological calibration, causal compensation or assay substitutability. The earlier frozen primary criterion remains failed.

## Consequential wording limits first

1. The four opposite-sign targets stable across maps occupy three sites, not four independent sites. The two c.39 substitutions share site52. Maps reuse measurements; these are not independent replications.
2. Every one of the212 site-omission calculations keeps its original predictions. Zero passes means that removing any single site's assessment contribution cannot make these fixed aggregate comparisons pass. It says nothing universal about refitting after deleting training data, other models, or biological absence of information.
3. The protocol and summary retain historical 'temporarily stopped' and 'pending separate review' language. This review supplies the final computational status without altering those historical bytes. They should not be presented as current unresolved numerical findings.
4. Well-conditioned observed fits and resolved empirical-range signs do not authenticate source preparation labels or convert low-bin enrichment into physical protein abundance. Both observations remain compatible with shared normalization error and proxy-specific effects.

## Contract and independent computation

Inputs are the three already authorized, pinned actual_v3 contrasts, predictions and evaluation artifacts. The reviewer did not import or execute analyze.py or any fitting/normalization function. No source counts, simulations or other outcomes were accessed. The script independently recomputes errors from frozen predictions and observed contrast centers using binary stored-float Fraction arithmetic, averages within sites and then equally across sites, checks all sensitivity subsets, and recomputes every fixed-prediction omission. The six diagnostic CSVs and full summary were checked. Every source and diagnostic input was hashed before and after.

Command from workspace root:

```
work/env/Scripts/python.exe -X utf8 work/research_cycle8/comt_followup/independent_review/verify.py --actual work/research_cycle8/comt_results/actual_v3 --followup work/research_cycle8/comt_followup/results_v1 --out NEW_REVIEW_OUTPUT
```

All four map/load cases retain125 targets,53 sites,113 SNVs, both regions and all five models. Their exact site-balanced aggregate and per-site scores match the archived results. The reviewer also independently reduces the frozen preparation log observations to medians and contrast endpoints. The largest difference for that alternative floating evaluation order is1.7763568394002505e-15, within the prespecified2e-12 arithmetic identity tolerance. All displayed CSV score, omission gain and range-radius comparisons have zero difference; exact fractions and flags use no tolerance.

Finite coverage:250 variant-stability rows,80 aggregate score rows,1,060 site/model error rows,212 omission rows,40 linear-fit diagnostic rows and1,500 extrapolation-flag rows. The flag CSV repeats each target across three relevant models; it is not1,500 independent cases. All training and assessment IDs are disjoint and their union is the complete target set. No site crosses a fold boundary.

## Omission and sign results

Code heavy, code mid, GEO heavy and GEO mid each have zero out of53 omission summaries passing all original comparisons. Joint-versus-ROI site wins/losses are25/28,24/29,25/28 and23/30 respectively; no exact ties. Omission gains versus ROI range respectively from -0.03239058546810816 to0.007027011371082626,0.008281668328910485 to0.04855760444326867,-0.05483369108621973 to -0.01005917922080258, and -0.005454210591788823 to0.037082475272365044. All848 baseline-specific omission gain entries were reconstructed, not only the ROI comparison.

Five code-map targets and four GEO-map targets have resolved opposite RNA/low-bin-proxy contrast signs. Their intersection is:

| Target | Site | Exclusive assay-log widening radius, both maps |
|---|---:|---:|
| c.33C>G |50|0.04823470901528504|
| c.39T>A |52|0.009140397875488038|
| c.39T>G |52|0.0016186334586221207|
| c.72G>C |63|0.07052958783613184|

Every other target remains in the table with its appropriate sign status and zero disagreement radius. This is empirical-range stability, not a confidence level or an effect-size equivalence test.

## Exact sensitivity algebra and counterexamples

For intervals A=[a-,a+] and B=[b-,b+], their Cartesian difference is[a--b+,a+-b-]. Widening each assay-log interval by delta changes these endpoints to a--b+-2delta and a+-b-+2delta. A contrast with strictly positive lower endpoint m remains strictly positive exactly when delta<m/2; the negative case uses m=-upper. For two opposite resolved signs take the smaller positive margin. At equality at least one contrast envelope includes zero and loses its strict resolved-sign label. This is an exact statement about the declared rectangular range enlargement. With extra joint constraints it remains sufficient but need not be sharp. No empirical covariance estimate or independent sampling assumption is required for this deterministic sufficient bound. Conversely, the absence of an authenticated error radius prevents a calibrated reliability, equivalence or noise-exclusion claim.

Known-answer boundary checks used rational margins1/8,1 and3/2. At delta=m/2 their lower margin is exactly0; at delta=m/2-1/1000 it is positive, and at delta=m/2+1/1000 it is negative. Common additive shifts of all raw assay logs cancel the contrasts, so a shared offset of arbitrary size can leave signs unchanged even though an unrestricted box of the same marginal widths loses resolution. A biological interpretation must preserve this distinction.

Opposite contrast signs do not identify compensation: take arbitrary assay-log proxies r=-1,d=0,f=-1. Then R=r-d=-1 and T=d-f=1, but no translation or protein degradation rate has been specified or measured. It is also possible for joint predictive gain to come entirely from cancelling shared RNA measurement error, as in the already sealed pre-outcome counterexample. These alternatives are not eliminated by the sign table or by site-heldout prediction.

Fixed-prediction omission is not refitting: with two training observations(0,0),(1,1), the fitted least-squares line predicts0 at0. Removing the first observation makes the slope unidentified without a separately specified fitting convention. Merely removing its assessment error leaves every original prediction unchanged. Thus the current omission summary cannot claim training-deletion robustness.

## Conditioning, rank and flags

Exact site weights are1/(training variants at that site). The reviewer forms weighted centered RNA/load Gram entries A,B,D using exact Fractions, checks determinant AD-B^2 and exact rank, and computes the two-feature condition number independently as lambda_max/sqrt(AD-B^2), where lambda_max=(A+D+sqrt((A-D)^2+4B^2))/2. This is the condition of the centered weighted predictor matrix, not the raw uncentered design including intercept. RNA-only uses one nonconstant feature and condition1. All40 exact Gram/rank sets match, including all stored coefficients and degeneracy labels.

Conditions range from1 to3.427336259029487. The Gram-eigen calculation differs from frozen condition numbers by at most4.440892098500626e-16 and from the independent-SVD CSV by at most2.220446049250313e-15; the recorded SVD-versus-stored maximum1.7763568394002505e-15 is reproduced. Numerical condition comparisons use relative/absolute1e-10. This tolerates numerical evaluation order, not biological measurement uncertainty. No folds were removed or refitted.

All1,500 flag rows match direct assessment-versus-training ranges, including strict outside comparisons and boundary ties. RNA-only models have no observed-load extrapolation diagnostic. RNA extrapolation counts are2 per case for each relevant model. Joint load extrapolation counts are2 in heavy and3 in mid under each map. These flags describe training-range position, not calibrated predictive uncertainty or proof of invalidity.

## Review execution history and limits

The independent harness initially encountered three schema assumptions before completing: the subset is named SNV rather than SNV_only; assessment IDs must be derived from prediction fold labels; CSV RNA_slope uses uppercase RNA. Each was corrected only in the reviewer adapter. These were not production/result mismatches, no source result was adjusted, and no output directory was created by the failed attempts. The completed receipt preserves the exact final code and input hashes. No model, threshold, fit, target, normalization or scientific criterion changed.

This audit validates the finite saved diagnostic computations and their bounded deterministic algebra. It does not reproduce RF training, source extraction, source metadata authentication, historical access order or independent biological preparation measurements. No new model was fitted, and the original failed hurdle was not rescued.
