# Frozen CNOT3 conditional prediction protocol, version 1

Status: scientifically frozen after independent synthetic implementation review. This protocol and the exact code/input manifest must be committed before the first target-outcome decoding. Phase-specific hash authorizations and the runtime guard are required. No gene-level half-life values were examined when these rules were fixed. Public summaries were already known, as recorded below. Commit dates are privacy-normalized and are not evidence of calendar-time registration; the archived phase records and parent-commit order document the declared sequence of work.

## Question and estimand

Does deposited CNOT3 IP/input log2 enrichment improve prediction of the deposited CNOT3-depletion log half-life beyond baseline log half-life, a fixed nonlinear baseline adjustment, complete codon composition and available reference/measurement covariates?

The primary quantity is the relative reduction in mean absolute error (MAE) on a fixed gene-group holdout, relative to the specifically named algorithms below. This is not a test of nonparametric conditional information. This is within-study prediction of processed gene-level measurements, not a causal pathway contribution, independent experimental replication, biological-replicate inference or an RNA structure mechanism. Reanalysis and original data are one source experiment. The reference is a transparent GENCODE48 proxy, not an authenticated expressed isoform.

The outcome is natural log(HL_CNOT3KO). Natural log(HL_WT) enters the predictors. Predicting the log fold change with baseline subtracted from both outcome and prediction has identical errors and is not independent corroboration. Conditioning avoids placing the noisy baseline explicitly on both sides of a change-score regression; it does not remove baseline measurement error, shared normalization, growth/dilution effects or confounding.

## Exposed information and reproduction gate

Published aggregate correlations, directions and sample counts are already known. Neither the prediction team nor root has read gene-level half-life values. The published-control reproduction will decode outcomes throughout its fixed cohort before modeling; subsequent holdout results will therefore be called locked, model-held-out prediction after aggregate reproduction, never wholly unexamined validation. Neither reproduction summaries nor holdout results may change the model hierarchy, partition, thresholds or exclusions.

The separately frozen reproduction reconstructs Figure 2 using archived STAT_IP, top1000 versus the remainder. Before fitting, require exactly5481 distinct eligible genes, positive finite paired half-lives, stored-ratio arithmetic agreement, positive finite SE_IP, STAT_IP arithmetic agreement and no exact tie spanning the top1000 boundary. Require median(high archived STAT_IP group HL_FC_LOG2) minus median(remainder HL_FC_LOG2)>0, the direction authenticated in the published figure. No published P-value threshold is a gate because caption and panel thresholds differ. Ratio and statistic arithmetic use abs(error)<=1e-10+1e-8*max(abs(stored),abs(computed)). Any failure stops the extension for source diagnosis; no alternate subset search is allowed. Published count/direction agreement is a reproduction, not new evidence. The original read filter cannot be rechecked from these columns and must remain a limitation.

Figure3 complete-case correlations and the PTMD overlap are secondary harmonized diagnostics because their exact source universes are not published. Discrepancies in these secondary quantities must be reported and cannot be repaired by changing cohorts. They do not by themselves authenticate or invalidate the explicitly defined Figure2 reconstruction. This distinction is fixed before values are opened.

In the reproduction discrepancy list, only PTMD_overlap_not_available_on_full_IP_cohort and published_PTMD_cutoff_tie_rule_unreported are secondary warnings. Every other code is a principal gate failure, including unrecognized codes. The existing negative-SE check conservatively covers every finite-IP candidate, even one later lacking usable half-lives; the zero-SE check covers the principal cohort. This broader conservative scope is retained. Count agreement cannot recover the unreported original read-coverage/missingness rules.

## Eligible universe and exclusions

Start with the frozen Figure2 reproduction cohort. Retain only exact stable gene IDs also present in the original SLAM conversion workbook after stripping a terminal numeric ENSG version, without aliases or duplicate resolution. First take the already fixed longest-valid-CDS representative from reference_v1. Exclude its gene if that selected transcript is not protein_coding, lacks complete CDS, or has mRNA_start_NF or mRNA_end_NF. Never select a shorter replacement. Require finite enrichment and positive finite SE_IP. Beyond the stated positivity/finiteness eligibility, no significance, expression, effect-size, half-life-range or outcome-magnitude/sign trimming is permitted. Exclude no gene merely because depletion reduces its estimated half-life.

Zero annotated UTR length is retained as zero. Its GC is imputed using the training-set median among nonempty UTRs, with a separate zero-length indicator. If an inner or outer training set has no nonempty UTR of the relevant kind, that fit fails; no arbitrary median is supplied. The annotated sequence may not be the expressed UTR. Other nonfinite reference features fail eligibility. Malformed reproduction rows and duplicate GeneIDs are integrity failures, not selectively removed observations. Report every exclusion and sequential counts. Reference coverage is eligible final genes divided by Figure2-cohort genes also in the original SLAM identifier inventory. Require at least80% coverage,2000 training genes,300 holdout genes,100 training groups and30 holdout groups. Every inner fold must have at least100 assessment genes and50 fitting-side training groups. Failure closes the primary test as insufficient design coverage, without a replacement split or more permissive cohort.

Let P be the principal operational cohort, S its intersection with the original stable-ID inventory, and A the final reference-eligible set. Report |P|,|S|,|A|,|S|/|P| and |A|/|S|; any required zero denominator fails closed. The80% rule applies to |A|/|S|, not to all original genes. Missingness is not assumed random. Grouping applies only to A and does not bridge related proteins outside A.

## Related-sequence grouping and fixed split

Use the selected complete proteins only to prevent near-identical representatives crossing partitions. Join proteins when full-length unit-cost Levenshtein distance d satisfies10*d <= max(lengths), and take transitive connected components. Exact identities are included. Use the proved lossless overlapping5-mer multiset bound, then exact global edit-distance scores. Broader paralogy, shared domains, regulatory dependence and experimental dependence are not eliminated by this rule.

Before real grouping, cap inverted-index pair updates at20000000, exact candidate pairs at50000, the sum of pairwise length products at3000000000 and elapsed grouping time at300seconds. The posting bound counts all shared5-mer pair encounters, including length-ineligible pairs. Elapsed time is checked during enumeration and after each native alignment call. An in-process clock cannot preempt a native call; no result completing after the cap is accepted. Stop on any cap, retaining an incomplete receipt; never drop edges, split a large component or replace the rule. These caps concern resource feasibility, not statistical selection.

Group ID is the lexicographically smallest stable GeneID in the component. Compute SHA256 of UTF-8('cnot3-focus-v1|outer|'+group_ID), interpret all32 bytes as an unsigned big-endian integer. Holdout if integer modulo5 equals0; otherwise training. Inner5-fold index is SHA256('cnot3-focus-v1|inner|'+group_ID) modulo5, on training only. No repartitioning after counts or outcomes. All sequence grouping and fold assignments are saved before model fitting.

## Features and models

The baseline basis contains log(HL_WT), its centered square and centered cube, plus positive-part hinges at its training10th,50th and90th percentiles. Center log(HL_WT) by the current fitting subset mean before constructing its powers; this preserves the polynomial penalty under a common lifetime-unit change. Quantiles use NumPy method='linear'. Each cross-validation fit learns its mean and knots from that inner training subset only. Duplicate knots are retained; ridge handles dependence.

The core covariates are log(CDS length), CDS GC fraction, log1p(UTR5 length), UTR5 GC, log1p(UTR3 length), UTR3 GC, the two zero-length UTR indicators and log(SE_IP). SE_IP is only an enrichment-precision proxy, not a half-life standard error. There is no valid decay-quality proxy in the recovered file. Missing half-life uncertainty, read coverage, source fit filters, censoring and growth corrections remain material limitations.

Compare six fixed models on exactly the same genes:

1. No change: predict log(HL_KO)=log(HL_WT).
2. Baseline linear: ordinary least squares with intercept and training-centered/population-scaled log(HL_WT), using SVD rcond=1e-12; a constant baseline column becomes zero.
3. Baseline/core: ridge on the nonlinear baseline basis and the core covariates.
4. Amino-acid/core: model3 plus20 amino-acid fractions.
5. Codon/core: model3 plus61 sense-codon fractions. This is the primary composition comparator.
6. Association: model5 plus the training-residualized measured LOG2_FC_IP defined below. This is the primary proposed model.

All ridge columns are centered and divided by their training population standard deviation; a constant column becomes zero. Center the response and leave the intercept unpenalized. Minimize sum squared training residuals plus lambda times squared standardized coefficients. Use lambda in {0.1,1,10,100,1000}; choose by pooled out-of-fold MAE across the fixed five training folds. Numerical ties within1e-12 prefer the larger lambda. Refit once on all training genes. No test-set tuning, feature selection, outlier deletion or new models.

AA fractions are aggregates of codon fractions. CDS GC is largely determined by codon composition and the terminal stop. These nested covariates are intentional benchmark definitions, not independent biological evidence. The published PTMD score is a fixed codon-frequency contrast and is not appended as new information.

To avoid apparent gains caused solely by duplicating a ridge-penalized covariate, residualize every proposed enrichment basis column using predictors only. Standardize that column by the current fitting subset mean and population standard deviation; a constant column becomes zero. Regress it on an intercept and the complete standardized model5 design using SVD least squares with relative rank cutoff rcond=1e-12. Apply that fitted projection to assessment features. If the training residual L2 norm is at most1e-10*sqrt(n_training), set the added column identically zero in training and assessment. Otherwise center it and divide by its training population standard deviation. Save means, scales, coefficients, rank, singular values, residual norm and zero flag. Every inner fold repeats this procedure using only its fitting subset. The no-SE sensitivity uses its corresponding no-SE comparator design; the nonlinear-enrichment sensitivity residualizes each proposed enrichment hinge as well as enrichment itself against the same comparator design.

This safeguard removes exact-spanned duplicate-feature effects at the declared numerical rank and fixed penalty. It does not remove nonlinear dependence, measurement error, unmeasured confounding or the possible amplification of a small noisy residual. Hyperparameters remain selected by training folds only. Synthetic exact-duplicate, near-duplicate and common lifetime-unit rescaling tests are mandatory.

## Locked prediction and endpoints

Write all six model predictions, training tuning records, preprocessing parameters, fold IDs and input/code hashes before invoking the holdout evaluator. The training API must accept no holdout targets. Commit the prediction file hash before scoring the holdout. Source outcome access by the prior reproduction remains recorded.

Primary success requires model6 MAE at least5% lower than model5 and at least5% lower than every other fixed baseline, with all eligibility and integrity gates passed. If any baseline MAE is zero, the relative improvement criterion is undefined or impossible and primary success is false; report absolute errors. This threshold is a practical pilot criterion, not a significance level or guaranteed biological importance. Report the actual incremental error whether positive or negative. If the primary test fails, do not tune against these holdout genes or call an alternative metric a successful validation.

Report MAE, RMSE, median absolute error,90th-percentile absolute error, mean signed error, and Spearman correlation of predictions and measured log half-life. Report both equal-gene and equal-component MAE, but success is defined by equal-gene MAE only. Save every heldout prediction/error and all model comparisons. Report error by the four intervals (-infinity,q25],(q25,q50],(q50,q75],(q75,infinity) defined by training baseline quartiles; retain empty intervals, with no subgroup success claim.

All quantiles use NumPy method='linear'. Equal-component MAE is the unweighted mean of component-specific MAEs. Constant-target or constant-prediction Spearman is undefined and stored as null. Success uses the exact floating comparison model6_MAE<=0.95*each_positive_baseline_MAE, with no additional scientific tolerance. Independent numerical verification uses declared calculation tolerances without changing the success rule.

No confidence interval or p-value will be labeled biological uncertainty because independent culture-level outcome estimates and half-life uncertainty are unavailable. Report the full range of primary improvement against the codon comparator after deleting each holdout component in turn, without refitting, as finite-cohort influence sensitivity. Report whether the sign and the5% criterion against every primary baseline survive every deletion. This is neither a confidence interval nor independent replication.

## Fixed sensitivity and descriptive kinetic comparisons

Two sensitivity models are permitted and cannot replace the primary endpoint: model6 excluding log(SE_IP), paired with model5 excluding it; and model6 adding enrichment hinges at its training10th,50th,90th percentiles, paired with the unchanged model5. Use the same split, tuning grid and training-only transformations. Report every result. They test precision adjustment and a fixed nonlinearity, not a feature search.

The nonlinear-enrichment basis is constructed from original LOG2_FC_IP and its fitting-subset quantiles before residualization; each linear/hinge candidate is then standardized and separately projected onto the unchanged comparator design. It is not a hinge transform of an already residualized score.

For the entire eligible cohort and separately for the holdout, report baseline half-life, log2 lifetime change, and the signed dimensionless rate contrast R=ln(2)*h_ref*(1/HL_WT-1/HL_KO), with h_ref the training median HL_WT. R is invariant to a common unknown positive time-unit multiplier; do not present it as an absolute per-minute pathway rate. The printed published axes say minutes, while the file-to-axis calibration and fitting conventions remain unverified. Report median,10th/90th percentiles, signs and Spearman correlations with each other and measured enrichment. Retain all negative contrasts. No pathway-specific rate attribution, rate-model fitting or uncalibrated propagation of nonexistent standard errors.

The reviewed kinetic bounds with residual activity and compensation are a separate identification result. They can explain why equal observations permit different pathway contributions, but no numerical biological bound is asserted without justified nuisance ranges. Shared residual activity across genes adds joint constraints; marginal sharpness is not simultaneous attainability.

## Verification and stop

Before outcomes, verify synthetic known answers for preprocessing isolation, ridge solutions, tuning, split invariance, grouping boundaries, negative controls, positive incremental association controls and shared-baseline coupling. Independent review must challenge the code and protocol. After predictions, a fresh implementation must verify cohort joins, split leakage, all predictions and metric arithmetic. Fix implementation defects transparently while preserving failed artifacts; a change to a scientific rule after exposure creates an exploratory successor, not this frozen test.

The execution runtime is the existing CPython3.12.14, NumPy2.5.4 and Biopython1.88 environment. Record actual versions and refuse a mismatch for the frozen initial run. Grouping is an exact-score calculation within its resource acceptance caps; no numerical model setting changes after results.

Stop after one principal reproduction and one locked model comparison with the specified sensitivities. Report a negative result as failure of this predictor under these measurements and models. Distinguish hypothesis failure from unusable source fields, inadequate coverage or nonidentifiability. No broad new biological branch follows automatically.

The initial execution invokes guarded_phase.py for each scientific phase. It verifies the exact script hash and frozen runtime before starting the child process. It supplies the runtime refusal gate for the prediction implementation. Source and phase inputs have separate mandatory hash checks.
