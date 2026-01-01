# PERSIST prediction protocol, revision 2

Status: proposed for independent review and root authorization. This package fixes the analysis before real folding, sequence feature calculation, development target access or holdout target access. It does not authorize execution. The exact machine-readable cohort, model settings, source hashes and feature columns are in `protocol.json`.

## Question and estimand

Does one newly computed ViennaRNA mean unpaired probability (MUP) improve prediction of the deposited `degradation_coef_in_cell` numbers, beyond five fixed comparisons, on transfer from the 22 prefix-310 constructs to the 15 prefix-930 constructs? Predictions and errors remain in the deposited numeric units. The unresolved absolute rate or half-life calibration is not repaired by this analysis. A known response transformation differing by row or source invalidates this endpoint and stops the comparison.

The unverified deposited `AUP.Vienna.FULL` column is excluded. No cell from it is decoded, reproduced, compared with the new feature or used to choose settings. This explicit replacement responds to missing provenance before outcome access. It is not an additional feature search.

## Immutable membership and dependence

Use all 37 longamp-listed rows with protein SHA-256 `dd0e8ca4f9a41c1fdd610ac0e48d953ff55be9483e2067db5e713804931fc148` and barcode-normalized context SHA-256 `e682ac5099bf1b6f4beb22f1f9817c89022c0bd31b279975a872546eddf726e9`. The ordered identifiers, exact sequence hashes and split are fixed in `metadata_settings.json` and copied into the protocol. Each full deposited RNA has 928 positions, the CDS has 666 positions including its terminal stop, and the barcode has nine positions. There is no membership revision from a calculated feature or response.

The development partition consists of the 15 connected components of the CDS Hamming-distance-at-most-ten graph on the 22 development rows. There are 12 singleton groups, groups {310054T3, 310055T3} and {310144T3, 310146T3}, and the group {310168T3, 310169T3, 310170T3, 310172T3, 310173T3, 310174T3}. Both exact membership and grouping were verified from permitted metadata. Prefix allocation is a design-source transfer, not random allocation. Longamp is a published processing flag potentially influenced by earlier outcomes. One reporter, one context, one study, shared design ancestry and unresolved biological preparation identity limit independence.

## One structural feature

Use the entire deposited `RNA.sequence`, including barcode and constant sequence. The only normalization is uppercase DNA letter T to RNA letter U. Other unexpected characters stop execution. Do not trim, substitute a construct, append cap or poly(A), or use experimental constraints.

The frozen engine is the installed ViennaRNA Python package 2.7.2 with double-precision partition functions and explicit Turner2004 parameter loading at 37 degrees Celsius. The complete binding-exposed effective model configuration, base-pair type matrix, reverse types, alias table, package file hashes and exported Turner parameters are archived. These include one linear strand; dangles 2; GU pairs allowed; isolated pairs allowed; minimum loop size 3; unlimited base-pair span; no G-quadruplex calculation; `betaScale=1`; `sfact=1.07`; salt 1.021 M; and no auxiliary constraints. The otherwise unexposed pair-dissimilarity field retains the pinned package default. No comparative-sequence fold is called.

For each sequence, create one fold compound, calculate its MFE once solely as a numerical scaling input, call the supported `exp_params_rescale(mfe)` method exactly once, and calculate the partition function once. There is no manual scale setter, fallback, alternative setting or retry. Record the resulting scale and both energy values as numerical diagnostics only. They are not predictors.

For N positions and upper-triangular probabilities p(i,j), the sole structural feature is

`MUP = 1 - (2/N) sum_(1 <= i < j <= N) p(i,j)`.

Validate this against the mean of `1 - sum_(j != i) p(i,j)` using the same probability matrix. The installed interface uses an (N+1)-square matrix, with unused index zero and a zero diagonal/lower triangle. Probability and marginal range tolerance is 1e-10; formula agreement tolerance is 1e-12. Values are not clipped. Require finite values. Record warnings from each isolated worker's output streams. Any warning, engine failure or invalid required value stops the entire comparison with a saved failure record. There is no row deletion, imputation or second model. The parameter export and exact package binaries are hash checked before every real stage.

## Predictors and models

RNA alphabet order is A,C,G,U. All frequencies have explicit denominators. The 21 composition predictors are four mononucleotide frequencies over N positions, 16 overlapping dinucleotide frequencies over N-1 pairs, and barcode GC fraction. The 61 sense-codon frequencies use the 221 CDS sense codons, excluding the terminal stop. Internal stops are invalid. The primary vector has 83 entries: 21 composition, 61 codon, and one MUP value. These fixed columns may be redundant; positive ridge penalties handle collinearity. No other stability, expression, response-derived or historical AUP field enters any model.

The six models, all reported, are:

1. Development response median.
2. Development response mean.
3. Ridge on the 21 composition columns.
4. Ridge on the 82 composition and codon columns.
5. Nearest development CDS under nucleotide Hamming distance, averaging all responses tied at minimum distance.
6. Primary: ridge on all 83 columns, adding only the new MUP.

For each ridge model, select lambda from {0.1,1,10,100} by deterministic leave-one-development-group-out prediction. Minimize total absolute held-fold error divided by 22. In each fold, exact constant training columns are removed and the remaining columns are centered and divided by their training population standard deviation (ddof 0). Holdout or omitted-group values never determine these moments, feature removal or coefficients. The objective is sum of squared residuals plus lambda times the squared coefficient norm. The intercept is unpenalized. Refit with the selected lambda using all 22 development rows and their own preparation.

A tuning tie has loss at most the best loss plus `1e-12 * (max(development response) - min(development response))`. Choose the largest tied lambda. When the development response range is zero, use exact loss equality. This response-range rule corrects the original absolute tie threshold and is homogeneous under positive response rescaling. It is applied to the complete development responses, not the holdout responses. No nested claim about out-of-sample development error is made from the same cross-validation used for tuning.

## Endpoint, coverage and decision

The primary endpoint is MAE on all 15 held-out rows. Operational success requires primary MAE no more than 0.9 times each of the five baseline MAEs, with every baseline MAE strictly positive. A zero baseline error makes positive gain impossible, including a tie at zero. Report all five gains, marking a zero-denominator gain undefined. This is a fixed-cohort effect-size decision, not a significance test.

Report every held-out prediction and error. Secondary diagnostics are RMSE and Spearman rank agreement with average ranks for exact ties; rank agreement is undefined for a constant vector. Keep negative and zero response numbers. Require all 37 predictors, all 22 development responses and then all 15 validation responses to be available and finite. A failure stops the complete-cohort analysis, preserving the gap. Do not inspect holdout finiteness before predictions are frozen.

Report all 15 leave-one-holdout-row-out MAE gain records and their ranges. These dependent influence diagnostics are not confidence intervals, new validation sets or grounds for excluding a row. No confidence interval, p-value, biological noise ceiling or bootstrap over supposedly independent genes is produced.

## Stages and access controls

1. Current stage: inspect permitted metadata and installed settings; test only synthetic matrices, arbitrary toy features and a 12-position synthetic fold. Hash the combined source table as bytes without decoding response cells. Preserve previous exposure incidents. No real sequence feature has been computed by this implementation task.
2. Root review and freeze: commit the exact protocol, source and tests, then provide an execution freeze with `status=authorized_by_root`, only `features` authorized, and all required input hashes. The candidate manifest in this package deliberately authorizes no stage.
3. Feature stage: compute the one rowwise structural feature and fixed sequence frequencies for all 37 rows. Open no response cells. Save complete diagnostics and the immutable feature artifact. Root checks this artifact and includes its hash in a separate training authorization.
4. Training stage: verify the complete feature artifact and all source hashes, then selectively decode only the 22 development values from the one endpoint column. Fit all six fixed models. Save all 15 predictions, every CV grid result, fold preparation, coefficients, nearest-neighbor ties and configuration hashes. Do not decode held-out response cells or old AUP cells.
5. Prediction freeze: root commits the complete predictions and writes a receipt containing repository, commit, Git path and SHA-256. This is separate from the stage execution authorization.
6. Validation stage: verify the prediction file hash, committed Git blob hash, implementation hash and protocol identity before selectively decoding the 15 holdout values. Report the frozen evaluation once. Preserve failure records and do not overwrite output directories.

Reading an entire byte stream for hashing or locating selected identifiers is distinct from decoding unrelated numeric cells. The selective reader decodes the header, identifier column and only the endpoint cells belonging to the authorized split. Poisoned unselected cells are covered by synthetic tests. These controls are inspectable workflow safeguards, not access restrictions against a person deliberately bypassing the code. Root remains responsible for the actual commit and authorization record.

## Interpretation and failure policy

A positive result would support incremental prediction by this explicitly modeled transcript feature on this selected transfer cohort. It would not establish measured cellular structure, causal CDS effects, a synonymous barrier, a complete mutational route, generalization to another gene or calibrated half-life prediction. Barcode composition and codon controls do not eliminate all selection or barcode effects. A negative result is retained as a failed predictive hypothesis conditional on this design.

Stop on ambiguous sequence or endpoint identity, hash mismatch, unexpected alphabet, required missing or nonfinite input, known source-specific endpoint conversion, unsupported engine settings, a folding warning, leakage, degenerate grouping or failed synthetic checks. No alternative feature, favorable subset or model retuning follows a validation failure. The earlier provenance failure for the old AUP field remains unresolved and separate.

The official ViennaRNA documentation supports the ensemble probability interpretation and the rescaling interface: [Python API](https://viennarna.readthedocs.io/en/latest/api_python.html). The installed version and binary hashes, rather than a mutable documentation page alone, define execution. The exact nucleotide-alias and pairing-array extents were checked against the [version 2.7.2 model header](https://raw.githubusercontent.com/ViennaRNA/ViennaRNA/v2.7.2/src/ViennaRNA/model.h).
