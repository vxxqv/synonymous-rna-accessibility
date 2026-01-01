# FAS measurement-contract review before outcomes

## Charter

Goal: decide whether the fixed72-genotype FAS pilot can answer a useful, bounded measurement question before any count/coverage/outcome access. Scope/object/subject: measured native-exon synonymous routes / endpoint and null-model identification / the proposed relative-enrichment envelope, forced-order test and source preparation contract. This is a descriptive adversarial review with inherited context, separately authored rather than blinded. Canon: originating Methods/GEO metadata and pinned author code first; native sequence audit second; mathematical arguments and synthetic checks third. Existing source eligibility is working knowledge; biological calibration remains incomplete.

Questions: Q1 choose enrichment or PSI; Q2 identify normalization and invalid-state rules; Q3 delimit conditional uncertainty and source confounding; Q4 challenge monotone-null and accessibility claims; Q5 decide whether outcome opening is worthwhile. No real count, phenotype, outcome missingness, read file or validation table was opened.

## Canonical sources

- Baeza-Centurion et al.2019, [10.1016/j.cell.2018.12.010](https://doi.org/10.1016/j.cell.2018.12.010), GSE111316 metadata and [Scaling_Law](https://github.com/lehner-lab/Scaling_Law) commit `b87064c39a553b0759d86c2b428813c58dc17ae2`.
- Prior source report `../measured_synonymy_sources/report.md`, native sequence enumeration, `sample_contract.json` and GEO protocol metadata.
- Author enrichment source blob `87838a729c4304ce3e38df87daaefb30c249119e` and calibration source blob `c3ee712c600e926a4a653410f860ac49c10b08bd`. Calibration source SHA256 `d05b9a172ceab2dc085379622db227f7cfe05fb1a0fd4116d5a64406947e41cb`.
- Root draft `../fas_pilot/protocol_draft_v1.md`, snapshot hash recorded in `review_receipt.json`. Recommendations refer to that exact reviewed snapshot, not subsequent edits.

## Findings ledger

| ID | Verdict and evidence | Consequence |
|---|---|---|
| F1 | **Accept relative included-RNA enrichment.** For positive reference/input counts, `R[v,b,t]=(O[v,b]/O[w,b])*(I[w,t]/I[v,t])`. Exact cancellation is shown below and checked with synthetic rational counts. | The target is recovery of variant-bearing included-RNA amplicons relative to plasmid input and human WT, not direct PSI. |
| F2 | **PSI adds unnecessary calibration dependence.** Masked source shows each published per-replicate PSI equals a common positive constant times its ES divided by ancestral ES. Nine formulae have identical positive calibration constants and matched replica indices; literal values were not displayed. | Within each replicate, relative-to-human-WT PSI cancels the calibration exactly. Averaging PSI before reference normalization can change replica weights; do not assume its mean equals mean relative enrichment. |
| F3 | **The27 envelope is descriptive.** Nine output transfections cross three technical input PCRs from the same library. Author ES uses the common median normalized input frequency; there are not27 independent biological preparations. | Report all scenarios and their outer box, with no confidence or population-coverage interpretation. |
| F4 | **Systematic recovery confounding is not identified.** Included-RNA output counts can depend on total RNA abundance, inclusion, RT recovery and PCR recovery. Input DNA and WT normalization do not cancel arbitrary genotype-dependent factors. | No rejection can by itself identify splicing-specific intercodon interaction or a protein/evolutionary mechanism. A concrete false-mechanism example appears below. |
| F5 | **The forced-order null has a valid finite representation under its assumptions.** Strictly separated closed intervals force strict latent order; finite strict inequalities can be rescaled. | Eight full-rank design columns include an intercept; the latent coefficients remain nonidentifiable under an arbitrary response. Exact certificates are necessary for a certified numerical contradiction. |
| F6 | **Two initial draft claim checks needed correction and root applied them during review.** Weakly nondecreasing transforms need not preserve rejection, and a global null contradiction plus a barrier need not show the particular barrier requires intercodon interaction. | Required distinctions and correction status below. These do not block the revised descriptive pilot. |
| F7 | **No independent validation cohort is presently authenticated.**2016 is the same exon/research programme and can inform the2019 design; potential18-genotype overlap is not untouched graph validation. | Hold back replicates only as conditional same-library repeatability, not independent biological validation. |

## Endpoint and normalization recommendation

Use **relative RNA-library enrichment**, explicitly an included-RNA recovery proxy for exon inclusion. Do not call it relative PSI unless the measurement assumptions required to equate enrichment and inclusion probability are stated. The source's selected RT-PCR calibration examples support assay correspondence for those examples, not uniform absence of bias over the fixed72.

For a separate input t and output b, the full-library normalized enrichment is

`E[v,b,t] = (O[v,b]/sum_u O[u,b]) / (I[v,t]/sum_u I[u,t])`.

Dividing by `E[w,b,t]` cancels both full-library totals, giving F1 exactly. No total is needed for this new separate-TR relative score. Preserve source totals only for the **separate historical author reproduction**, which uses medians of input frequencies. In that reproduction the actual R code first aligns every file to TR1's sequence index and uses `na.rm` sums and medians. Therefore reproduce both that alignment and the resulting totals; all-file totals need not be identical if row universes differ. Record the discrepancy instead of silently replacing the author's total or omitting missing selected vertices. Do not normalize within72 and pretend to reproduce the author calculation.

There is also a useful enclosure argument: if all input frequencies are positive, the ratio `median_t a[t] / median_t b[t]` lies between `min_t a[t]/b[t]` and `max_t a[t]/b[t]`. Indeed, componentwise `a[t] <= M*b[t]` implies `median(a) <= M*median(b)`; the lower bound is analogous. Consequently the source's relative median-input score for each output lies within that output's three separate-input relative scores, provided the same sequence universe/count definitions and finite positive inputs apply. This proves an observed-estimator enclosure, not coverage of a true biological quantity. A synthetic example corroborates the algebra; the argument, rather than that example, establishes the general relation.

The masked PSI code also shows outcome-derived low-noise filtering and full-source mean/SD summaries. The new pilot must not borrow those filters or present the calibration data as training-only or withheld. No calibration outcome values were exposed in this review. All248 numeric tokens were masked before display, including numbers in comments/identifiers. Source code was never executed. The sealed byte copy is local provenance only and is not proposed public-package content.

### Invalid counts and zeros

The draft's fixed72/unknown policy is appropriate, with these executable interpretations:

1. Each scenario requires positive `O[w,b]`, `I[w,t]` and `I[v,t]`, and authenticated records. An input WT zero cannot be interpreted as a valid zero numerator. A missing or zero reference invalidates the affected complete scenario. Do not force WT=1 when its reference observation is undefined.
2. If any required scenario for v is undefined, the primary union-envelope for v is `[0,infinity]`, rather than the range after discarding that scenario. Keep every vertex and record the cause. This may make a null untestable; that is a legitimate information limitation.
3. A finite observed zero output with valid positive reference/input gives observed R=0 without a pseudocount. It does not establish zero true inclusion or zero biological recovery. Repeated observed zeros remain sampling-limited even if the empirical range is a singleton.
4. Negative, malformed or noninteger counts and conflicting duplicate rows stop the pilot. An omitted row is not zero. No coverage-based replacement or favorable subcube after access.

The minimum/maximum is a coordinatewise rectangular outer envelope of the27 scenario vectors. It permits combinations that no preparation produced and may therefore be conservative for this finite model-compatibility question. It is not a probabilistically calibrated identified set for the true phenotype. Report both the full scenario vectors and outer-box status; do not count scenarios as independent trials or use their range as a confidence interval.

## Statistical and biological identification limits

Write a conceptual included-RNA measurement as `O[v,b] proportional to D[v,b]*A[v,b]*p[v,b]*q[v,b]`, where D is the effective transfected template representation, A is transcript abundance per template, p is inclusion probability and q combines RT/PCR/sequence recovery. Sequenced input I may estimate source plasmid representation rather than the exact effective intracellular D. Normalization can cancel shared scales, but arbitrary genotype-dependent A, D/I or q remains. The same intended protein does not imply identical RNA stability, processing or assay recovery. Constant primer binding sites and equal amplicon length help control some technical variation but do not establish identical RT/amplicon recovery. Source-specific mechanisms that remain uncalibrated cannot be repaired by more computation.

**Known-answer counterexample.** On a two-codon binary square, let true inclusion be1/2 for every genotype. Let persistent recovery factors be1 at00/11 and1/5 at01/10, with identical source DNA. Observed scores are1/2,1/10,1/10,1/2. For `z=a+b*x+c*y`, a common nondecreasing response would require `c<0` from `y00>y01` and `c>0` from `y11>y10`. The observed score null is rejected even though true inclusion has no interaction. Repeating the same systematic bias across all nine preparations and all three inputs leaves a singleton envelope and does not resolve the ambiguity. The exact synthetic receipt retains this example.

Thus a negative codon-unary/order result can identify **incompatibility of the deposited included-RNA enrichment ordering** with that named model class. It cannot distinguish intercodon splicing regulation from genotype-dependent RNA abundance/recovery, an incorrect common-response assumption or other stable assay effects. These are alternative explanations, not evidence the data are bad. No causal contribution, endogenous isoform abundance, actual protein production, apoptosis, clinical or evolutionary claim follows.

Replicate variation is still informative for descriptive repeatability. Resampling a whole output genotype vector, rather than treating vertices/squares as independent, would preserve a preparation's shared noise if a later inferential protocol were justified. The current source does not authenticate independent culture-day or cloning-batch blocks, and no bootstrap confidence claim is recommended. Three technical input PCRs cannot create independent biological validation. The present use of all nine transfections for an envelope means none remains a validation arm for a later method.

## Mathematical contract and required draft corrections

For each v, let the nonempty closed interval be `[L[v],U[v]]`, with finite nonnegative L and possibly infinite U. If `U[i] < L[j]`, any compatible nondecreasing response must have `z[i] < z[j]`. Conversely, suppose every such strict order is respected by `z=X*b`. Sort vertices into equal-z blocks and assign each block the maximum L among that block and all preceding blocks. This value cannot exceed any U in its block: otherwise an earlier/same-block lower bound would force a strict reverse z order. These assignments define a common nondecreasing response inside every interval. Equal z values receive one value. This proves the claimed equivalence within these interval/response assumptions; it does not prove a biological model.

With finitely many strict latent orders, a feasible coefficient vector has a positive minimum margin across the forced pairs. Multiplying all coefficients by its reciprocal gives the margin-one inequalities. If no pair is forced, a constant latent vector and a compatible constant response exist because all intervals overlap pairwise and one-dimensional closed intervals have a common intersection. Infinite upper bounds contribute no constraint of the form `U[i] < L[j]`; they can still be the high vertex of an order forced by a finite U elsewhere.

For infeasibility, a nonnegative normalized dual witness with `A^T*lambda=0` and `sum(lambda)=1` contradicts `A*b>=1`. The exact rational verification should also check the forced pairs against exact rational count comparisons. Floating rounding near ties is not an effect-size criterion. A feasible witness needs response assignments checked against every interval, not merely a small LP residual. Unknown intervals do not become failed genotypes.

Required distinctions raised during review of `protocol_draft_v1.md`:

- **Transformation check, synthetic-check paragraph.** Replace “positive monotone transformation” with **“one common strictly increasing transformation, applied consistently to interval endpoints and all vertices”**. A weakly nondecreasing constant transform can erase all forced orders and convert rejection to acceptance. Positive output values alone do not imply strict monotonicity. Logs at zero need explicit domain handling, not an unannounced pseudocount.
- **Co-occurrence claim, opening question and accessibility interpretation.** State that global null rejection and a preparation-consistent barrier are two descriptive findings; their co-occurrence does not establish that the particular barrier requires an intercodon interaction. For example, on `P3 x Q2`, set `y(k,x)=q(k)+r(x)`, `q=(10,0,10)`, and `r00=r11=1`, `r01=r10=0`. The Q2 slice rejects the common-monotone codon-unary null. Yet the barrier between `(0,00)` and `(2,00)` already comes from the P3 codon gate. A separate targeted necessity argument would be required for an interaction-dependent barrier claim. Do not add an outcome-selected attribution test.
- **Design rank, null paragraph.** Prefer “an eight-column full-rank categorical design” to “eight identifiable columns.” Free h leaves latent coefficients nonidentifiable up to at least positive affine changes; it is not an eight-parameter model overall. Exact design rank and contrast dimension are computational checks, not biological identifiability.
- **Conditional certification wording.** Use “certified relative to the declared empirical envelope and72-state design” whenever a lower/upper certificate is called conservative or certified. No interval coverage for true mean inclusion has been established. “Preparation-consistent” is safer than unqualified “reproducible.”

**Correction verification.** Root applied the transformation, co-occurrence and coefficient-identifiability fixes while this review was being completed. They are present in the held final `reviewed_protocol_snapshot.md`, lines33,49,59, and are closed rather than outstanding blockers. The initial read preceded those edits; its byte hash was not captured and is not invented. The final snapshot hash is recorded in the receipt. The remaining conditional-certification terminology is an interpretation requirement already supported by the draft's explicit empirical-envelope/graph limitations, not a reason to add an inferential claim or another analysis.

The graph may be disconnected at a threshold under arbitrary codon-unary effects because each Leu factor is a three-state path. A nucleotide-additive failure is therefore not enough. The primary codon-unary/common-response null is the appropriate stronger comparator. The inspection establishes no new foundational result or novelty priority beyond the sources already reviewed.

## Scope filter

Retained a conditional, descriptive included-RNA enrichment question. Removed any inference that WT-normalized counts identify direct PSI, that27 combinations are independent replicates, or that a repeated assay bias is eliminated by an empirical envelope. No new outcome branch, calibration fit or predictive benchmark is proposed.

## Object filter

Relative enrichment eliminates unnecessary global total/calibration scales, but not variant-specific recovery. Separate median-source reproduction remains a historical computation and cannot silently replace the new scenario definition. Graph and protein identity remain the full fixed72 native reference class, not an outcome-selected cube.2016 remains selected same-exon complementary material with possible design influence.

## Conclusions

**Qualified GO for one frozen descriptive72-genotype pilot under the corrected snapshot and conditional claim boundary above.** F1-F7 support an inspectable, falsifiable model-compatibility and bounded-accessibility calculation. A feasible null, lack of certified obstruction or uninformative envelopes are useful possible results and cannot trigger retuning.

The opening is scientifically worthwhile only if the deliverable accepts the limited estimand: which orderings and route restrictions are supported by the deposited included-RNA enrichment observations across the declared normalization/preparation envelope, and whether that envelope excludes arbitrary codon-unary effects under a common nondecreasing response. It is **not** worthwhile as a proposed resolution of splicing mechanism or calibrated biological significance, because this source alone cannot distinguish the counterexample explanations above. Stop that stronger claim rather than describing a computation as a measurement advance.

Suggested result wording, conditional on a genuine positive computation: “No assignment within the predeclared coordinatewise empirical enrichment ranges was compatible with the named codon-unary/common-response model on this72-genotype native-reference translation class. This is conditional on the deposited assay and envelope, and does not identify a splicing-specific mechanism.” Report any accessibility result separately and qualify its graph/envelope scope.

## Next steps / Risks

Freeze the endpoint, full72 membership, all27 scenarios, original-normalization reproduction, zero/unknown rules, exact order comparisons, solver certificate verification and all-endpoint report. Independently review the implementation and synthetic cases before any observed counts. No separate new calibrated uncertainty model is recommended in this task.2016 outcomes, selected source calibration observations and held-back transfections cannot be called untouched independent validation without their own exposure/design audit.

### Exposure boundary

Existing source Methods, metadata, pinned enrichment code and this draft were read. The only newly retrieved artifact was the exact known author calibration `.R` source; all248 numerical tokens were masked before model display, including embedded calibration values. A separate source-property check disclosed only that all nine per-replicate formulae use the same positive coefficient and matched indices. No literal coefficient, calibration outcome value or count/PSI table was displayed. No author code was executed. The only executed calculations used synthetic rational counts and the stated synthetic confounding example. No observed count, coverage, outcome missingness or2016 validation material was read.
