# Independent reserved-validation verdict

**Computational verification passed with zero mismatches. Scientific validation did not pass.** The frozen independent adapter returned `verified_complete_arm`, meaning the selected export comparison completed; this status does not mean adequate scientific coverage or a successful validation result.

Only6 of24 selected reserved squares were eligible, representing6 of16 groups and12 diagonal comparisons. The unchanged minimum is10 eligible squares spanning8 groups. Reserved coverage fails, and therefore the frozen two-arm descriptive milestone is not met. Neither pooling development with reserved observations nor highlighting a surviving positive can repair that failure. The conditional prediction screen also fails on both coverage and numerical-performance criteria.

## Reconstructed result

| Quantity | Independently verified reserved result |
|---|---:|
| Eligible squares / selected squares |6 /24|
| Eligible groups / original groups |6 /16|
| Eligible diagonals |12|
| Point joint-only obstructions |2 /12, or16.6667%|
| All12 deposited-series-consistent diagonals |1|
| Fixed-low-end empirical-range certificates |1|
| Scorable conditional prediction cases / groups |6 /6|

The adapter independently checked12 diagonal descriptions,144 whole-series scenarios and108 exact critical-threshold cells, including point, range, missing-component and numerical-tolerance classifications. It reconstructed312 source-mean diagnostic records, retaining missing/non-comparable cases, with no source-mean provenance hold. There were12,049 exact leaf comparisons and1,761 numerical-identity comparisons. Scientific labels, coverage, thresholds and tolerance values match exactly; arithmetic identities agree within the frozen absolute/relative2e-12 comparison tolerance. No scientific threshold was adjusted.

The complete point-positive inventory is:

| Square / group | Point observed component | Component outcomes over12 series | Empirical ranges |
|---|---|---|---|
|6 /277|Complete; joint-only obstructed|8 joint-only;4 single-trait disconnected|Local and component inconclusive|
|20 /46|Complete; joint-only obstructed|12 joint-only|Local and component fixed-low-end certificates pass|

Square20 is a **category B finite processed-data observation**, conditional on the stated empirical intervals and fixed low-end thresholds. Its certificate does not calibrate those intervals, imply universal endpoint-relative obstruction over every perturbation, establish an independently replicated biological effect, or rule out unmeasured synonymous detours. The certificate does not override inadequate arm coverage.

## Prediction failure

All three frozen methods and both traits remain reported. The table displays independently reconstructed equal-group errors with rounded precision; exact stored values are in the receipt.

| Method | Protein MAE | Protein RMSE | Growth MAE | Growth RMSE |
|---|---:|---:|---:|---:|
|Additive completion|4.90774965|10.03345194|0.01817211|0.02913943|
|Adjacent-corner mean|4.50915931|5.94545291|0.01602234|0.02174055|
|Three-corner mean|7.30261847|10.32357211|0.01987772|0.02642655|

Relative to adjacent-corner mean, additive MAE is worse by8.8396% for protein and13.4173% for growth. Relative to three-corner mean, its improvements are32.7947% for protein and8.5805% for growth. Thus three of the four ten-percent comparisons fail. There is no success obtained by choosing the favorable trait or baseline. Prediction coverage also fails the original10-case/eight-group minimum.

Two positive and four negative prediction labels occur. Additive completion has pooled confusion TP2/FN0/FP0/TN4 and group-balanced classification error0; each baseline has TP0/FN2/FP0/TN4 and error1/3. Additive balanced accuracy1 on these six cases does not override either failed screen condition and supplies no calibrated future discrimination claim.

## Scope, counterexample and conclusion

The overall frozen validation question is **category D: insufficient coverage**, with separately retained category B local observations and a failed predictive benchmark. This is not a category E source failure, because the declared source gate passed. It is not a biological null, a practically-small-effect exclusion or category A independent replication. Shared source processing, reporter/library, growth founder/reference and uncertain deposited-column generating revision remain unremoved dependencies.

The strongest mechanistic counterexample is exactly additive antagonism: on Q2, F=(x1-x2,-x1+x2), endpoints00/11 are(0,0) and intermediates have opposing signs. Each trait admits a route while their intersection admits none, although both interaction contrasts are zero. The same graph observation therefore cannot establish epistasis, RNA folding, natural selection or biological function improvement. Common strictly increasing transformations of final trait summaries preserve the primary classification while altering additive contrasts and projected rank. A sampled component's missing path also cannot exclude a route through unmeasured genotypes.

The positive finite examples may be described transparently, including native margins, every series result, interval assumptions and measured graph scope. They cannot be presented as successful reserved validation or as resolution of biological uncertainty. No frozen model, endpoint, coverage minimum or threshold should be changed to rescue this result. Supplemental amendment diagnostics retain their separate verification requirement.

## Post-outcome analytic limitation of the classification baselines

This limitation was recognized after the frozen result comparison. It does not change the frozen baseline, labels, numerical thresholds or pass/fail rules, and is not a new theorem-priority claim.

**Elementary coordinatewise-convexity lemma, exact arithmetic.** On exactly Q2, let s,t be the designated opposite endpoints, z the known intermediate and c the predicted withheld intermediate. For each trait j set tau_j=min(s_j,t_j), using weak acceptance. Suppose

```text
c_j = lambda_s,j s_j + lambda_t,j t_j + lambda_z,j z_j,
lambda_s,j,lambda_t,j,lambda_z,j >= 0,
lambda_s,j+lambda_t,j+lambda_z,j = 1.
```

The weights may differ among traits; common weights are a special case. If z_j>=tau_j then every input in this convex combination is at least tau_j, so c_j>=tau_j. Thus z's set of acceptable traits is contained in c's set. If every trait individually connects s to t, then each trait passes at z or c. The containment forces it to pass at c in either case. Consequently c passes all traits and the route s,c,t is jointly acceptable. **A convex predictor of this form cannot produce joint-only obstruction on this square.** Equality and zero weights are included. If the individual-connectivity premise fails, the conclusion is only that the prediction is not joint-only, not that it must be jointly connected. No larger-graph or unrestricted sequence-space characterization is asserted.

Adjacent-corner mean has lambda_z=0 and both endpoint weights1/2, so c always passes every trait and explicitly provides a joint route. Three-corner mean has all weights1/3. It need not itself pass every trait, but its accepted-trait set still contains that of z, so it also cannot give joint-only obstruction in exact arithmetic.

The minimal exact witness distinguishes additive completion's negative coefficient. Let s=t=(0,0), z=(1,-1). Three-corner mean c=(1/3,-1/3) fails trait2; both intermediates fail that trait, so it predicts individual disconnection. Adjacent-corner mean c=(0,0) gives a joint route. Additive completion c=s+t-z=(-1,1) has a negative coefficient on z; its intermediate pass patterns oppose those of z and produce a joint-only obstruction. This is compatible with an exactly additive trait field and does not establish epistasis.

There is a finite-precision qualification: the exact lemma must not be advertised as a universal guarantee for every floating implementation at weak ties. For synthetic float endpoints s=(0.7,0),t=(0.7,1),z=(0.7,-0.1), ordinary three-term addition and division yield c_1=0.6999999999999998<0.7. The tiny rounding error can break acceptance-set containment even though the exact rational mean ties at0.7. This is a synthetic arithmetic warning, not a detected mismatch in the authorized results. Their baseline labels were independently verified to be all-negative in both arms. The exact example and this numerical distinction are recorded in `reserved_run/convex_baseline_analytic_check.json`.

Accordingly the observed classification advantage over these baselines is structurally weak: in the intended exact geometry, the comparison is against two all-negative joint-only classifiers. The additive model's perfect finite confusion table remains reported, but it is not independent evidence of mechanism or strong predictive discrimination. The continuous-error benchmark remains separate and failed as already specified. No replacement comparator is fitted or retrospectively added.

Stop this bounded independent review here.

## Exact verification provenance

Authorization SHA256: `2adfab21315150672e35c19c5991b223914222c3a4070f6d462e5504897991d8`.

Frozen executing adapter SHA256: `6b16d291a8c94804403c1fdea6d281d035892aa8d2d2eeeacd53d81da20d050e`.

Authorized evaluation SHA256: `0f36e4d762b1e80e87859e700f4da2655c6ecf99c29e0560d9d6dcd9546a5521`.

The unchanged adapter was executed from the repository with the exact supplied authorization and wrote only the new `work/joint_function_cycle5/observed_result_review/reserved_run` directory. All ten input hashes matched before and after. The receipt is `reserved_run/receipt.json`; its full hash list is reproduced below. This verifies selected-export consistency and independent arithmetic, not raw source generation, historical access order or biological calibration. No raw archive, other outcome source, external service, frozen code or repository file was accessed or changed beyond the explicitly authorized reads.

Receipt SHA256: `b5d9886b79d49eebdfe652b65435f0b3193d25577af467bc14d26493fee68dd9`.

| Tested input | SHA256, unchanged before/after |
|---|---|
|evaluation|`0f36e4d762b1e80e87859e700f4da2655c6ecf99c29e0560d9d6dcd9546a5521`|
|unsealed_inputs|`50fdb6099cdbffc91802753888cbe66a9e7d46e83aff20e2a5c83260ce74d0c7`|
|predictions|`716cf3fde5eb5d7e2c9361bee06880fc752d1511df2f9e377ea720b1a77feac8`|
|prediction_seal|`e12a78e8051d9d5c2d6590c6a0b62f68bee3fb9003b3cb4e72019f7ed9719f35`|
|known_inputs|`30cb26844153cb080dfb7b11056b3a0682c647a731608bd688ce623d8fc53e86`|
|split|`394a0c1aa23785c290b31c147d199f406cf46b269daa9159a98c419398b05a4e`|
|components|`90772dc57ab6dfdadf737e72cc0be11783245299d1ec0e9e589663a64ee1dadd`|
|edges|`5d7d5f3b1160d7a8338da3f39369720fdc92c6105b47efbd902eeede48a30edd`|
|adapter|`6b16d291a8c94804403c1fdea6d281d035892aa8d2d2eeeacd53d81da20d050e`|
|reference|`64a5e9040a77d81dc779583d1d3ea5894a5aebd5718dbf3466281803f48ce375`|
