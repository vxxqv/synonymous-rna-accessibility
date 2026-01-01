# COMT A1 independent pre-outcome methods review

**Verdict: PASS for the held source/estimand amendment and synthetic model-kernel building blocks, subject to one criterion clarification. HOLD for any observed count decoding or outcome execution.** The selective reader, integrated normalizer, preparation joins, staged outer-fold driver and guards are explicitly unfinished and were not certified by this review. No source phenotype, count body, confidence limit or source simulated-training value was opened.

This is a separately authored review with inherited project context, not a blinded external referee. Reviewed held candidates are identified in `review_manifest.json`. Historical `protocol_draft.md` was read to understand the abandoned design; its blocked state and sequential PAVA-residual comparator are not assessed as the final method.

## Required correction before freeze

**Protocol v2, line 57:** make the success quantifier unambiguous: require the 10% MAE reduction against **all four baselines under each of the two documented FACS grouping maps**, if that is the intended rule. The present sentence explicitly quantifies coverage over maps, but “a favorable result that reverses” could allow a positive 9% gain in the second map after a 10% gain in the first. A direction reversal and failure to meet a fixed magnitude threshold are different. This is a prospective wording/driver requirement, not an observed result or permission to tune a threshold.

## Surviving scientific question and source limits

The full construct/reference and exact nucleotide joins support a descriptive held-out-codon-site prediction question for the deposited Flag-COMT-moxGFP system. The 125 source-ascertained synonymous identifiers at 53 sites, with a fixed 113-SNV sensitivity, are defensible fixed ascertainment. They already inherit the source's filtering and four-assay identifier availability. They are not a new random sample or an outcome-independent source cohort. Site folds keep all alternative codons together, protecting the fitted calibration from seeing another allele at an assessment codon. They do not hold out biological preparations, regions, cell lines or source preprocessing.

The target T=d-f correctly reverses low-bin enrichment. It means reduced representation in the source's low-moxGFP/mCherry gate, not positive protein amount or a calibrated mean fluorescence. RNA/polysome and cytometry induction durations differ. Neither transcript annotation nor full-fusion reference synonymy authenticates all mature isoforms actually produced. These limits are preserved in the held protocol and source contract.

Four stable cell lines per ROI are source-supported, but the 20 conflicting FACS sample labels and repeated preparations are not resolved by this pilot. Reporting both named maps is appropriate sensitivity, not proof that either is correct or that they bracket every possible truth. Medians of technical measurements followed by four biological-group summaries are permissible defined descriptive statistics. Missing required groups must remain incomplete; incomplete groups cannot silently become extra technical evidence or a three-line design.

The Cartesian preparation envelopes are descriptive unpaired bounds. They need no assumed cross-assay pairing, but have no calibrated coverage, cannot supply simultaneous biological confidence, and may admit jointly impossible R/T/L combinations because r and d recur. Resolved sign contrasts are sufficient statements about the declared ranges, not probabilities of molecular direction. Existing ALDEx2 marginal intervals do not repair this gap.

## Shared-reference counterexample and claim boundary

For R=r-d, T=d-f, L=h-r, independent centered raw noise gives Cov(R,T)=-Var(d), Cov(R,L)=-Var(r). This algebra creates association without a layer-specific mechanism. The saved exact synthetic control is stronger: true RNA=heavy=u, true heavy/RNA load=0, measured RNA=u+e, d=0, f=-u. Then R=u+e, L=-e, T=u. Across u in {-1,1}, e in {-2,2}, the constrained joint model has a=b=1 and zero MAE, versus nondecreasing PAVA MAE 1/2. A large load-comparator gain can therefore arise from correcting RNA measurement noise while true load is constant.

The v2 protocol now includes this control and explicitly excludes causal compensation, translation flux and independent biological-information claims. A positive result would support only predictive gain for the named derived proxy and preprocessing algorithm on held-out sites. A null would concern those named algorithms and ascertainment, not impossibility of predicting protein from RNA. Neither alone establishes biological assay equivalence; no source-supported practical substitution margin is available.

## Model capability and mathematical checks

The five fixed models are substantially more defensible than the discarded sequential residual model. The joint constrained RNA-linear plus unrestricted load-linear fit contains the constrained RNA-linear baseline. The ROI-specific weighted-median baseline guards against simple region-level offsets; pooled median remains separately reported. Nondecreasing weighted PAVA is a capable **named** RNA-only benchmark, not every RNA-only predictor and not the MAE-optimal monotone fit: it minimizes training squared error and evaluates MAE. Its clamping and interpolation differ from linear extrapolation, as declared.

For centered sufficient statistics xx, zz, xz, xy, zy, the determinant xx*zz-xz^2 is nonnegative. At positive determinant the stated normal equations give the unconstrained optimum. If its RNA slope is negative, the boundary a=0, b=zy/zz satisfies the constrained KKT conditions: the one-sided RNA residual derivative has the required sign. If load varies but is collinear with RNA, the prescribed a=0 load fit spans the identifiable predictor direction and is feasible. Constant-load and both-constant branches reduce correctly to the constrained univariate fit or weighted response mean. The predetermined collinear convention affects extrapolation but is not outcome-selected. Rank/condition and out-of-training-range diagnostics must remain visible.

I independently compared PAVA with exhaustive contiguous-partition minimum weighted SSE, consolidating exact predictor ties, and verified constrained least-squares KKT conditions, site weights, ROI medians and threshold boundaries: **166 checks passed**. The copied, byte-unchanged author harness also passed **550 synthetic checks** in this review's new directory. These tests corroborate the code; the convex projection/KKT arguments justify the general optimum characterization. No author numeric function served as the independent reference oracle.

Site-balanced scoring correctly averages within site then across sites. The future driver must aggregate all assessment-site errors directly, rather than average five fold errors equally when folds have different site counts. Shared model coverage and all dropouts must be preserved. Report each ROI and the fixed SNV sensitivity regardless of direction; neither a successful subset nor the mid-load sensitivity may rescue a failed primary.

## New runtime and source-normalization amendment

The amendment explicitly defines a new fixed-runtime reconstruction, preserving the failed historical-equivalence gate. It pins R/packages, helpers, seed, RF tuning/defaults, source truth/factor conventions, rejection-ID thresholds and row order. It retains the executable ROI2 exact-half rejection rule rather than silently substituting the prose majority description. It preserves raw REF/ALT before source type classification, and requires background filtering before W reconstruction. These are necessary distinctions: full-codon rewriting can change the source classifier type, and a denominator using only the 125 targets is a different statistic.

I verified the held runtime manifest against **3,677 task-library file hashes and the recorded Rscript hash**. No mismatch was found. This verifies the manifested local bytes, not historical RF equivalence, independent official archive signatures, every system/OS binary or behavior on real source calls. The author RF smoke tests were inspected as receipts, not rerun or relabelled as independent real-data verification. Source helper bodies are pinned and allowlisted; the future caller must enforce RF-only use and fresh controlled startup.

## Remaining mandatory integration gates

- Run RF/reader stages in a fresh `Rscript --vanilla` process with the frozen local library, runtime/source hashes and sequential settings checked before source-row decoding. The helper environment inherits `globalenv()` (`source_helpers.R:29`), so uncontrolled startup/global bindings must not become a hidden algorithm input.
- Finish and independently test source member/schema/field allowlists, raw representation to canonical target joins, full background allele lineage, per-position depth uniqueness, MNP collapse and negative-control subtraction. Invalid/negative W, conflicting identities, incompatible features or nonfinite arithmetic must stop. Do not use a favorable alternative filter state after evaluation.
- Preserve complete source lineage, both grouping maps and biological-group completeness. Report RNA/target-only eligibility separately from the common RNA/load/target model cohort so loss due to the added covariate is explicit. Do not infer missingness at this stage; no values have been opened.
- Freeze outer folds, five model fits, all-map success quantifiers, no-optimization runtime guards, overwrite refusal, coefficient/rank/range diagnostics and complete per-variant/site/fold scoring. Keep source reconstruction sealed before prediction fitting and evaluation. Only prediction calibration is training-local; measurement normalization and source-ascertainment history are not an independently trained preprocessing pipeline.

These are already acknowledged unfinished components, not newly discovered numerical failures. Until they pass, this review's PASS cannot be used as authorization to run observed data. The whole-pipeline status remains HOLD.

## Evidence and limitations

`shared_reference_checks.json` records exact synthetic algebra; `independent_kernel_checks.json` records the 166 references and runtime-file checks; `held_kernel_replay/proxy_synthetic_receipt.json` preserves the author harness replay. `source_and_estimand_findings.md` contains the earlier preflight observations. The manifest records every reviewed candidate and test hash. No new primary retrieval, source outcome/count read, model selection on data, external write, commit or publication occurred. No unresolved source pairing or biological calibration was claimed resolved.
