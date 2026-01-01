# Independent CNOT3 pre-outcome protocol review

## Verdict

**REVISE the reviewed draft before freeze. The bounded estimand is defensible after the required clarifications and accepted amendments below.** It can compare fixed algorithms for prediction of deposited gene-level point estimates in this study. It cannot establish conditional information beyond all functions of the controls, a direct CNOT3 rate contribution, biological-replicate uncertainty, transcript-specific composition control, independent experimental validation or a mechanism.

This is a separately authored review with inherited context, not a blinded external review. No real half-life, association, standard-error or trajectory outcome cell was opened. Only the draft, source/reproduction/annotation records and source code were inspected. The official published Figure2 was re-viewed for its median direction. All new calculations use explicitly synthetic arrays. The root revised the shared draft during review. The snapshot/hash in `reviewed_sources.json` captures the amended80-line version, which was read again before this report closed. Scientific amendments described below are present there; they are not a certification of an implementation that has not yet been inspected.

## Required changes before freeze

### 1. Strengthen the baseline and remove the exact-span ridge artifact

The initial draft required the following amendments. They are now present in snapshot lines37,45-56,70 and must be preserved in the frozen implementation:

- Add square and cube of **training-centered** log(HL_WT) to the existing linear/hinge basis for models3-6 and their sensitivities. Learn centering, knots, imputation and scaling only inside the applicable inner-training or outer-training subset. Center before taking powers. Raw polynomial powers change the ridge penalty under a common conversion of lifetime units.
- Replace raw enrichment in primary model6 by its training-only least-squares residual against intercept plus the complete standardized model5 design. Fit the projection separately inside every inner fold and again on outer training; apply it without refitting to the matching assessment/holdout genes. Standardize enrichment before projection; use an explicit SVD relative cutoff, accepted as1e-12. If the training residual Euclidean norm is at most1e-10*sqrt(n_training), make the added column identically zero in both training and assessment. Otherwise center/scale its residual using training parameters. Save rank, singular values, coefficient, residual norm and all transformations. Do not amplify roundoff residuals to unit variance.
- The no-SE sensitivity must remove SE from both its comparator and its enrichment residualization design. The nonlinear-enrichment sensitivity must residualize each proposed candidate column against the same matching comparator, rather than reintroducing raw spanned columns. Its candidate-column and projection order must be explicit, including whether the raw linear enrichment column is projected before constructing hinges. Recommended: construct linear/hinge columns from original enrichment using training knots, standardize each, and project each onto the unchanged model5 design.

These changes require no additional outcome-driven feature choice or tuning branch. They do not guarantee nonparametric adjustment. Residualization removes the component in the declared numerical column space, not arbitrary nonlinear dependence or measurement confounding. Small genuine residuals can still be dominated by noise. Below-cutoff singular directions are not part of the retained numerical space; avoid claiming an exact symbolic span theorem for an SVD implementation.

Why this is consequential: with synthetic T=log baseline, X=T^2 and Y=T^2, X contains no information not already determined by T. The draft's three-hinge comparator nevertheless has holdout MAE0.1613562, while adding X gives0.0002551, a99.84% apparent gain. The accepted centered polynomial comparator gives0.0002557. These are known-answer model-class examples, not biological findings or a full valid-cohort fixture.

The second example uses identical standardized columns. Splitting one coefficient over k identical copies changes its effective penalty from lambda*b^2 to lambda*b^2/k. In the synthetic check, adding a fourth copy to three copies decreases holdout MAE24.99%, with no new feature information. The complete primary success gate includes other baselines, but that does not remove this general penalty failure mode. The accepted residualization makes an exactly spanned candidate zero and addresses it directly. For an exactly orthogonal residual, the ridge normal equations separate at fixed lambda; a tuned model can still select a different lambda, so report this as an algorithm comparison.

In a separate synthetic unit-conversion check, raw polynomial powers change predictions by up to0.73184 on the log scale, whereas training-centered powers differ by at most7.55e-15. The SVD example's exact-span residual norm is7.66e-14 versus zero threshold4.90e-9 and must be suppressed. Full receipts and independently authored code are retained next to this review.

### 2. Freeze operational reproduction, not a claim of recovered original eligibility

Draft lines17-25 are implementable with the separately supplied reproduction protocol, but must cite its exact frozen hash and make the scopes unambiguous.

The inspected `reproduce.py` first requires unique nonempty IDs and finite archived STAT_IP, LOG2_FC_IP and SE_IP, then finite positive paired HL and finite stored HL_FC_LOG2. It retains zero/negative SE and arithmetic mismatches as gate issues rather than selecting them away to recover5481. This is the right distinction. Count agreement alone cannot authenticate the absent >=10-read filter or the authors' unpublished missingness rule. The estimator is restricted to this operational complete-case reconstruction.

Both arithmetic tests already use atol1e-10 and rtol1e-8, with tolerance atol+rtol*max(abs(expected),abs(stored)). Copy or explicitly reference these exact rules. Use log2(KO)-log2(WT), not a potentially overflowing ratio. Archived STAT_IP remains the ranking variable. The source describes STAT as fold-change divided by SE; consistency with LOG2_FC_IP/SE_IP is expected only when those exported fields represent the same coefficient and scale. Disagreement is a source-identity/processing gate, not evidence of no association. Do not recompute the ranking, drop inconsistent rows or infer a shrinkage convention from a matching result.

The code currently places primary and secondary warnings into one `discrepancies` list and always holds for root review. The final gate must distinguish principal integrity failures from secondary PTMD/correlation-universe warnings, as promised by draft line19. Otherwise a missing PTMD score or PTMD cutoff tie can silently stop a valid principal comparison. Also specify the scope of negative-SE checks: `ip_candidates` flags any finite-IP candidate with negative SE even when that row later lacks usable HL, whereas zero-SE is flagged only in the final cohort. A broader conservative stop is permissible if declared; it must not be mistaken for a cohort-only rule.

Make the direction gate executable. Source prose says relative lifetime increment is bigger in the high-IP group. The official Figure2B shows a higher pink top1000 median, and its caption identifies the box's internal line as median. Thus require median(high stored HL_FC_LOG2)-median(rest)>0. This is a figure-supported direction check, not a directional KS test. Do not require an exact published P-value: caption says P<0.0001 while the panel says p<0.0005. Preserve this source discrepancy. No workbook outcomes were consulted for this recommendation.

For denominator reporting use explicit sets: P=operational principal cohort, expected|P|=5481; S=P intersect original SLAM stable-ID inventory after the one permitted version normalization; A=S passing the predeclared pinned-reference and finite-feature rules. Report |P|,|S|,|A|,|S|/|P| and|A|/|S|, every sequential exclusion and whether each ID would have been train/holdout. The80% rule is|A|/|S| and does not describe coverage of all original genes or prove missing-at-random selection. Zero denominators fail closed. Positive/finiteness eligibility is an outcome-validity condition; say no additional magnitude/sign trimming rather than imply literally no outcome-dependent exclusions.

### 3. Remove representative and imputation ambiguity

Draft line23 should say: first take the already fixed longest-valid-CDS representative from the reference table; exclude its gene if that selected transcript is not protein_coding or has either mRNA incompleteness flag; never choose a shorter substitute. The current phrase 'longest eligible NMD/LoF representative' confuses selection with the downstream restriction. Reference v1 explicitly permits NMD/LoF during initial longest-valid selection, then records the selected biotype. Preserve this two-stage rule and its hashes.

Zero annotated UTR length is not a missing gene and is not evidence of no expressed UTR. Retain it as declared, impute GC training-only and preserve the indicator. If an inner/outer training subset has no nonempty UTR values from which to calculate a median, stop with a preprocessing feasibility failure; do not invent a global or holdout-derived fallback. Save imputation medians, constant-column flags, knots, centering, scaling and the exact CDS/codon/stop conventions. No folding feature is included. There is no authenticated measured-folding variable, and this is neither an RNA-structure adjustment nor a test of its absence.

### 4. Define zero/error and missing-metric behavior

The initial draft needed an explicit positive-denominator rule; amended snapshot line62 now includes it. A relative gain is undefined if a comparator MAE is zero; require strictly positive comparator MAE and model6 MAE<=0.95*comparator MAE for every required comparator. A zero-error comparator makes primary success unavailable, not an automatic pass from0<=0. Preserve raw losses and the undefined ratio. Apply the same rule after each deletion and for sensitivities. State numerical comparison tolerances once, without tuning them after observed differences.

For baseline OLS use a fixed least-squares/rank convention, including a constant baseline. Constant-target/prediction Spearman is undefined, not zero. Empty baseline-quartile groups are reported empty. Pin quantile definitions for all90th/10th-percentile summaries and quartile interval boundary/tie rules; they may be repeated or empty with tied baselines. Define equal-component MAE as the arithmetic mean of within-component MAEs, not pooled gene errors or size-weighted component means. Influence deletion changes evaluation membership only and never refits or retunes; it is not a confidence interval. No calibrated uncertainty is recoverable from gene-level SE_IP.

## Six comparator and sensitivity review

| Definition | Verdict and necessary scope |
| --- | --- |
|1, no change, logKO prediction=logWT|Valid fixed reference. It does not assume true control/KO equality and is not fitted.|
|2, intercept+logWT OLS|Valid affine calibration comparator, using training only; fix the rank-degenerate convention.|
|3, baseline/core ridge|Valid algorithm baseline after centered polynomial amendment. Available IP precision and reference UTRs are proxies, not full measurement-quality adjustment.|
|4, model3+20 AA fractions|Valid fixed representation benchmark. With ridge it is not merely an unpenalized nested-space test.|
|5, model3+61 sense-codon fractions|Valid primary comparator for a pinned representative-CDS feature set. The effective penalty differs from model4; a performance difference does not isolate synonymous causality.|
|6, model5+IP residual|Accept the replacement described above. It is still a fixed-algorithm predictive comparison. Raw enrichment must not remain an undeclared alternative primary.|
|No-SE paired sensitivity|Valid if both comparator and residualizer drop SE, use the identical genes/split and refit only the predeclared training rules. It cannot resolve unmeasured decay-quality or abundance confounding.|
|Nonlinear-IP paired sensitivity|Valid if linear and hinge candidate columns are defined and residualized training-only, with the same comparator, tuning grid and reporting. It is an additional functional-form challenge, not independent validation and cannot rescue primary failure.|

All ridge predictions require the exact sum-squared-residual objective, not mean loss with the same numerical lambda. Pooled inner-fold MAE must weight assessment genes equally, using one prediction per gene; it is different from mean of fold MAEs when fold sizes differ. Standardization, polynomial centering, imputation and enrichment projection all belong inside each fold. The declared absolute tie tolerance1e-12 is acceptable in the fixed log-loss units, with larger lambda selected; save every candidate loss and chosen lambda.

## Grouping, resources and exposure

The near-identical-protein rule, transitive components, fixed SHA split and failure-on-cap policy are coherent. The supplied5-mer proof uses multiplicities, including the short-protein/identity exception, and makes the candidate bound necessary rather than sufficient. Every surviving candidate still needs exact global edit distance. No outcome or resource-driven deletion, edge dropping, splitting or retry with a weaker threshold is allowed.

Specify that grouping is applied to A, the final frozen eligible set. It protects direct passing pairs inside A and their transitive closure, not all external paralog paths. Fix source/protein normalization and lexical ordering, count pair updates on unique protein records as proposed, and define elapsed timer start/end and software/runtime before running. The300s cap is a platform-dependent feasibility gate, not a portable complexity guarantee. The count, length-product and group/fold minima may fail. Such failure is an unexecuted/underpowered-by-design test, never a biological null. No real grouping was rerun for this review.

The draft correctly records prior aggregate exposure and that published reproduction exposes all candidate cohort half-lives before fitting. Hence 'locked model-held-out after reproduction' is accurate; 'untouched', 'prospective', 'independent validation' or an unexamined holdout is not. Freeze all scientific rules and code before reproduction, then preserve membership/projection/prediction hashes before holdout scoring. Root review after reproduction can decide stop or a transparently labeled successor; it cannot silently repair scientific rules and keep the same confirmatory label.

## Rate scale and inference boundary

For positive B,K and training h_ref>0, R=ln(2)*h_ref*(1/B-1/K) is dimensionless. Replacing both lifetime columns by c times themselves also multiplies h_ref by c, so R is invariant for one common c>0. Its sign agrees with log(K/B); its rank need not. A gene-specific or condition-specific rescaling, undocumented censoring or a different growth correction invalidates this simple invariance. This is not a calibrated pathway-specific per-minute rate. The actual published axes say minutes; the unresolved issue is file-to-axis calibration and estimator conventions, not absence of a unit label.

Natural-log KO prediction with baseline as predictor avoids the explicit algebraic shared-baseline response but does not remove noisy-baseline proxy effects. An association correlated with latent baseline can improve prediction conditional on an error-prone measured baseline without predicting a distinct perturbation mechanism. Missing half-life errors, growth correction, true expressed isoforms, source coverage and assay abundance remain unmeasured. SE_IP is not decay uncertainty. No numerical nuisance-bound result or biological confidence interval follows from this protocol.

Acceptable positive wording: 'The frozen association-residual algorithm reduced held-out prediction error relative to the prespecified comparator algorithms for processed half-life estimates within this experiment and reference proxy.' A negative result is failure of that algorithm/threshold, not proof that CNOT3 association, structure or the mechanism is irrelevant. No priority claim follows from this protocol review.

## Optional improvements, not freeze blockers

Report effective ranks, enrichment-residual training variance and fraction suppressed as zero for each fold. Report train/holdout baseline and source-quality predictor distributions without using them to alter the split. These diagnose extrapolation and weak residual signal; they do not create new success criteria. A larger flexible learner, alternative grouping database, sequence-feature screen or raw trajectory refit is unnecessary for this bounded question and should not be added in response to results.

## Verification and remaining gates

Synthetic counterexamples were independently implemented without importing the author model/reproduction implementation. They support the specific algebraic risks and revisions above; they do not test unseen production code. Source selection code was inspected, not executed on protected data. No failed real-data analysis or alternate cohort search occurred.

Before authorization: integrate these rules, hash the exact protocol and implementation, add synthetic tests for centered-unit invariance, duplicate/near-duplicate candidate suppression, training-only residualization/imputation, constant/empty cases, positivity/zero-loss rules and all primary-versus-secondary gate scopes. Then obtain the requested independent implementation review. After that, one operational reproduction and one locked model comparison are defensible within the stated limits. Stop on unresolved principal source or resource failures rather than interpret them as hypothesis outcomes.

## Amendment verification at closure

The captured draft SHA256 is066f30406e0054b99c9248a3d6d5c9404f2ec794bb7313c328c4d7b6f7ae29db. Read-back confirmed the centered polynomial baseline (line37), primary residualization and numerical zero gate (lines48,54-56), matched no-SE/nonlinear sensitivity projection (lines54,70), positive-denominator failure (line62), empty-imputation failure (line25), explicit quartile intervals (line64), cap timing semantics (line31), and figure-supported median-direction/tolerance gate (line17). These resolve the principal original scientific-design objections. Remaining pre-freeze requirements are the exact cohort/primary-versus-secondary gate and SE-diagnostic scope, representative-selection wording, complete metric/OLS conventions and independent implementation tests. The amended design is acceptable in principle for the explicitly named within-study algorithm comparison, conditional on those executable requirements. No new scientific branch or real-outcome analysis is recommended.
