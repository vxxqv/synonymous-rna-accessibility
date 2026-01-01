# Independent review of PERSIST endpoint proposal

**Verdict: restrict.** Prediction of the deposited `degradation_coef_in_cell` numbers is scientifically defensible as a retrospective test of a published RNA-persistence proxy. Absolute half-life calibration is not necessary for that narrow numerical prediction question. The proposed split can test transfer between two deposited design-source strata within one reporter/context. It does not establish a biological mechanism, independent-study generalization, or an independently replicated decay-rate relationship.

**The proposal is not ready for outcome evaluation as written.** The absolute tie tolerance breaks its asserted scale-equivariance, and exact provenance of `AUP.Vienna.FULL` is not authenticated by the source material recovered here. The proposal itself makes feature authentication a prerequisite. Retain that gate; do not replace an explicit requirement by an inference from the column name.

## Endpoint interpretation

The author R source fits time coefficients from a shared Long-amplicon cellular measurement matrix and writes their negatives. The construct table has the corresponding named cellular coefficient field. The general Methods identify a cellular time-course RNA endpoint; they do not identify a genotype fitness or protein-output endpoint. Predicting a deposited response field in declared numeric units is a valid, inspectable artifact-level estimand despite uncertain absolute calibration.

Call it the “deposited cellular degradation coefficient” or “deposited cellular RNA-decay proxy”. Do not call its values natural-log decay constants, accurately calibrated hours-based half-lives, or true rates inferred up to a proven common scale. Historical transfer/renaming into the combined construct table is not reproduced. The visible log2-processing/ln(2)-conversion inconsistency is a reason to withhold absolute physical interpretation, not a mathematical reason to prohibit prediction of the stored numbers.

The source suggests common processing: one selected sample matrix, one design and one `lmFit` call across selected constructs. There is no inspected evidence of different response units by source. That observation does **not** establish that all measurement distortions are one common affine transform of true biological rates. Spike-in normalization, finite-count transformations, regression weighting, selection and preparation differences can have genotype-dependent effects. Positive-affine invariance only protects a hypothetical common change of numerical scale. A discovered source-specific transformation or a column assembled from noncomparable endpoints would invalidate this pooled numerical transfer question and must stop it.

## Affine-invariance argument and required local correction

For a fixed feature matrix, a ridge model minimizing

`sum_i (y_i-b-X_i beta)^2 + lambda ||beta||^2`

with unpenalized intercept and fixed lambda obeys `y' = a*y+c`, `beta'=a*beta`, `b'=a*b+c` for a>0. The entire objective is multiplied by a^2. Training-only feature scaling and constant-feature dropping are unchanged by a response transformation. Hence its predictions transform affinely. Development mean, median and the average of nearest-neighbour ties have the same property. Absolute errors and RMSE scale by a, while MAE ratios and ranks are unchanged when defined.

Selecting lambda by the minimum development MAE also has this property **if its tie rule is homogeneous**. The proposed absolute `1e-12` near-tie rule is not. If two candidate losses differ by `5e-13`, the rule can regard them as tied and prefer the larger penalty; multiplying all responses by ten changes their difference to `5e-12`, which can change the selected penalty. Such an arbitrarily small difference can always be obtained by rescaling any nonzero candidate-loss difference. The global invariance assertion therefore fails for the stated selection procedure, even though it will matter numerically only near a tie.

Use exact loss ordering followed by a fixed deterministic rule on exact ties, or a tolerance proportional to a homogeneous, translation-invariant training response scale, for example `1e-12*(max(y_development)-min(y_development))`. Declare the zero-range convention. A fixed positive floor in response units would reintroduce the problem. Floating implementations should promise verified numerical agreement over a declared scale range, not exact equivariance for arbitrary underflow/overflow conditions.

Also define the zero-error case. “At least10 percent smaller MAE” must not be implemented only as `MAE_primary <= 0.9*MAE_baseline`: if both errors are zero, that inequality passes despite no improvement. If any comparator has zero MAE, strict relative improvement over every comparator is impossible. Report it as such and leave the ratio undefined where appropriate. Undefined rank correlation for constant vectors is a diagnostic, not a reason to change the cohort.

## Design-source split, dependence and leakage

The exact protein/context/flag membership is already established by sequence metadata. The22 prefix-310 development rows and15 prefix-930 holdout rows have no duplicate CDS across the split; minimum cross-source CDS Hamming distance is118. The declared threshold10 development graph is nondegenerate:15 connected groups, consisting of twelve singletons, two two-row groups and one six-row group. Leaving one such group out leaves at least16 training rows. These facts were recomputed using metadata only and are recorded in `grouping_metadata.json`.

Those groups control literal near-neighbour overlap under the chosen metric. They do not authenticate independent design seeds or designers. A crowd submission can derive from an algorithm, and several algorithmic designs can share a seed. “Different deposited design-source strata” is defensible; “independent design methods”, “independent families”, “independent biological validation” or “independent study” is not established. The threshold and direction are prespecified, so do not search them after outcomes.

The source strata were chosen by sequence-design procedures. The primary paper describes AUP as a design objective for some designs. Consequently the held-out covariate distribution is selected, not a random sample of possible synonymous sequences. This does not leak the new cellular response labels into the proposed regression, but it can restrict range, change correlations and make transfer conditional on the published design process. The endpoint remains worth testing as a modest transfer challenge, including a negative result.

Exact barcode-normalized context matching removes large nominal-context differences; it does not make actual barcodes identical. Barcodes contribute to full-RNA composition and can affect transcript behavior, processing or structure. Full-RNA mono/dinucleotide frequencies plus barcode GC are useful nuisance controls, but they cannot isolate a causal coding-region effect from barcode effects. An incremental full-transcript AUP prediction result may reflect CDS or barcode-dependent structure, and other source-linked differences remain possible. Keep it a prediction of the actually assayed reporter constructs.

The deposited `longamp` flag is also potentially outcome-informed: in the earlier identifier audit its192 TRUE IDs exactly matched the192 deposited solution-output IDs, while the author source describes filtering anomalous solution time courses. The proposed37-row estimand is therefore conditional on deposited selection. Preserving all37 intended IDs and stopping on missing required values is appropriate. It does not remove selection bias or restore an outcome-independent sampling frame.

Training-only standardization, unpenalized intercepts, fixed features, positive penalties and grouped development tuning are reasonable. Codon frequencies are necessarily linearly dependent; positive ridge penalties with an intercept still make coefficients estimable. Do not add an unnecessary full-column-rank exclusion. Tune each declared ridge model using the same frozen procedure, without comparing alternative hyperparameter ranges on the15 held-out labels. Development CV selects penalties; its reported error is not an unbiased performance estimate after that selection. The15-row transfer comparison is the primary evaluation.

## Feature provenance: authenticated meaning versus unresolved field mapping

The primary [Wayment-Steele et al. AUP paper](https://doi.org/10.1093/nar/gkab764), equations4 and6, defines unpaired probability from a predicted base-pairing probability matrix, then averages it across the sequence. Thus AUP is a sequence-model-derived ensemble summary, not a measured structural value or, by definition, a response-trained score.

The primary [PERSIST paper](https://doi.org/10.1038/s41467-022-28776-w), Methods “Design algorithms/analysis”, says ensemble calculations used LinearFold/LinearPartition with ViennaRNA, CONTRAfold and EternaFold parameters. This matters: “Vienna” need not mean an exact conventional ViennaRNA partition function. Its construct-table header contains `AUP.Vienna.FULL`, but the pinned repository tree contains no producer for that column. The inspected R source produces cellular/ribosome endpoints, not the structural predictors.

The following narrow authentication status is supported:

| Item | Status |
| --- | --- |
| Generic AUP definition and sequence-based conceptual provenance | Authenticated from primary definition |
| Generic use of Vienna-parameter ensemble calculations in PERSIST | Authenticated from primary Methods |
| Existence and spelling of deposited column | Authenticated from header only |
| Exact producer mapping that column to full `RNA.sequence`, including experimental barcode/flanks | Unverified |
| Meaning of FULL versus PCR at this exact field, modeled poly(A)/cap treatment and included positions | Unverified |
| Historical parameter version, LinearPartition settings/beam, temperature and package convention | Unverified |
| Evidence that this exact deposited field was postprocessed using cellular outcomes | None observed; absence is not authenticated producer provenance |

AUP's generic equation does not establish all of those field-specific choices. The exact-name web search did not locate a producer. A read-only GitHub connector search for the generic helper name in `DasLab/arnie` returned no results; that is not evidence that no relevant source exists. The author-linked software page returned403, and a public supplementary descriptor request failed DNS resolution before any HTTP response. Neither route was retried or bypassed. The descriptor is unavailable through that attempted route, not demonstrated absent. No supplementary outcome workbook was opened.

The exact feature gate therefore remains closed. A producer script or a primary legend explicitly mapping the field to input sequence/model/options could close it. No new structure calculation or predictor substitution was authorized by this review. If provenance remains unresolved, the current proposal's specified primary model must stop. A different explicitly opaque “deposited feature-column prediction” question would require a separate root decision and freeze; it would not authenticate the structural interpretation.

## Executable protocol requirements before outcomes

1. Record the immutable37 IDs and22/15 allocation, protein/context hashes, the15 development groups, exact source files/hashes, feature list, model conventions, loss and penalty grid. Correct the affine-tie rule and the zero-baseline-error rule before freezing.
2. Authenticate the exact AUP field as above, or stop this primary comparison. Freeze its input/provenance and characterize it as predicted full-transcript structure, not measured structure. Confirm that no protected cellular response was used to generate or select the feature in this protocol.
3. Separate loading stages. First load the22 development targets and permitted features only; do not validate/inspect the15 held-out target values during development. Fix penalties, all preprocessing and all15 predictions, then write their hashes. Only afterward load the15 targets and apply the complete-cohort finite-value gate. A missing/invalid held-out target stops the defined comparison; it does not trigger a smaller analysis.
4. Implement synthetic leakage tests: changing, permuting or replacing held-out targets must not alter any selected penalty, preprocessing or prediction; fold preparation must be independent of its held-fold target and of holdout feature moments. Test positive-affine equivariance of complete tuning and prediction over stated scales, including negative responses, zero-range development targets and tied losses. These are proposed checks, not executed here.
5. Evaluate all six declared models once, preserve each error, and apply the frozen10-percent operational criterion without recentering on holdout labels. Report per-row influence diagnostics as dependent sensitivity calculations, not confidence intervals. Retain a failure without changing source direction, cohort, grouping distance, endpoint, transformation or model list.

Unresolved biological-repeat identity, the selected sampling frame and one-study/one-protein scope remain limitations even after these gates pass. They prohibit a biological noise ceiling or broad generalization claim; they do not inherently prohibit this transparent deposited-response transfer test.

## Sources, exposure and independence

Pinned local inputs and hashes are in `source_manifest.json`. Primary source URLs are the PERSIST DOI above, its [public Europe PMC XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8940940/fullTextXML), the [author repository](https://github.com/barnalab/persist), and the AUP definition paper above. Hashes identify the locally recovered source files; no hash of an unrecovered page is invented.

No cellular response value, solution-table row, ribosome-load value or model fit was opened or evaluated in this review. The first local grep inadvertently returned an already-published Methods paragraph with selected24-reporter aggregate correlations; a later AUP-definition extraction returned aggregate correlation text and design-objective descriptions. Web search also returned published aggregate prediction claims from related papers. Those source-text exposures are recorded; they were not used to choose the endpoint/split or to fit a model. The earlier solution first-row incident remains in the prior audit and is not relabelled untouched.

This is a separate review of the root proposal, with inherited project context and access to the previous audit, not a blinded review. All writes are confined to `work/persist_cycle3/proposal_review/`. No root proposal, source file, repository, manuscript, freeze, account or publication was modified. No fitting, cellular outcome evaluation, author contact, commit, or subagent launch occurred. The review stops at this restricted verdict and explicit feature-provenance gate.
