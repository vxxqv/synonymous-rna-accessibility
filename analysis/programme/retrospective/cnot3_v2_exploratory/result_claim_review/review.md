# Independent CNOT3 result and claim review

**Verdict: pass after one kinetic-quantifier correction, now verified in the current manuscript and report. No outstanding numerical or substantive CNOT3 claim correction was found.** This is a separately authored, source-aware review with inherited context. It is not a blinded external review or an additional biological validation.

## Quantitative trace

All109 recorded checks pass in `checks.json`. Headline quantities were traced to the exact archived evaluation, independent verification, reproduction/source-diagnosis summaries, cohort/design records and committed-prediction receipt. Numeric fields in the model, prediction, influence and kinetic CSV exports match their archived counterparts exactly. This covers nine models,10,008 model/gene error rows,1,109 component deletions and5,416 kinetic rows. The arithmetic CSV separately retains exactly the49 original failed identifiers, and every row fails the old tolerance while passing the stated ten-digit diagnostic. Its raw-source/Decimal columns were not independently regenerated in this claim audit. The complete Markdown reports in `outputs/CNOT3_Results_and_Decisions.md` and computational `result_report.md` are byte-identical, SHA256 `408b10ec3d464124e89d37db4cb68d656ac88fc333cc5ee36cdb68fa7c7dfeec`; all27 model-table entries match their declared rounding.

| Headline | Evidence and verdict |
| --- | --- |
|Exploratory operational population5480, published5481 unresolved|Independent source summaries and frozen reconstruction agree:22,848 deposited rows,7,359 finite-IP candidates,1,879 missing all three HL fields,5,480 retained. The precision diagnostic does not establish historical preprocessing. Pass.|
|5,416 reference-eligible genes;4,304/1,112 train/holdout;4,280/1,109 groups|Eligibility/design records agree, with98.8321168% reference coverage. Independent cohort audit reports all72 candidate distances checked. Pass.|
|Primary MAE0.3393159972, codon0.3397002408, AA0.3394583794|Relative gains recomputed as0.1131125498% and0.0419439468%. AA is the strongest primary MAE baseline. All displayed rounded values agree. Pass.|
|Frozen5% against every baseline fails|Archived `primary_success=false`; gains also fail against baseline/core, AA and codon. Passing against no-change/linear comparators cannot substitute. Abstract, main section, supplement and report retain this failure. Pass.|
|No-SE/nonlinear sensitivities|Matched gains0.0995179486% and0.3193327682%; they are fixed sensitivities and do not replace primary failure. Pass.|
|Primary RMSE worse than codon|0.4551274724 versus0.4550418171. Pass.|
|Deletion gain range|0.0782871314% to0.1397746958% against codon; positive codon gain survives every deletion, positive gain over every baseline does not. Influence interpretation is accurate. Pass.|
|Independent numerical agreement|Receipt passes59,268 checks/240,144 values, maximum3.3972824554e-13, rounded upward to3.4e-13 in text. Its evaluation hash matches the reviewed evaluation. Pass.|
|Kinetic description|All5,416 rows retained;912 negative contrasts. Rate/fold-change rho0.826797, baseline correlations-0.581621/-0.141675 agree. Derived endpoints share measurements; they are not independent evidence. Pass.|

The reviewed prediction receipt identifies unscored operational predictions, and the root-verified commit receipt binds the same prediction SHA256 to commit `b709472c473e1f7ff53a20d8c7fe28ec97b43bf1`. This review checks the archived receipt and consistency, not historical Git access order. It performs no Git operation.

## Required correction found and resolved

The initial `cnot3_extension.tex` sentence near line143 said that unrestricted background shift alone makes every p in[0,k0] compatible and then chose r=0, immediately after specifying potentially restrictive bounds r in[rL,rU]. With a physical perturbed background, that statement needs r=0 to remain admissible.

Counterexample: k0=2,k1=1 and rL=rU=1 imply k1=(b+delta)+p with b+delta>=0, hence p<=1. Candidate p=2 is impossible for every delta despite p being in[0,k0]. This does **not** falsify the displayed sharp-set inequalities, which correctly include rL*p<=k1.

Root corrected the manuscript to: 'If r=0 remains admissible and delta is unrestricted'. The user-facing/computational results report has the same qualifier and explicitly warns against arbitrary positive lower residual-activity bounds. Read-back confirmed the correction. The abstract now also states that the0.113% gain fails the5% criterion, resolving the suggested clarity improvement.

## Kinetic proof and scope

The displayed identified set is correct for fixed comparable rates, nonnegative b,p, independent Cartesian nuisance ranges and the stated optional physical background. At fixed p, w=(1-r)p must lie in both[(1-rU)p,(1-rL)p] and[D+dL,D+dU]. Their intersection gives the two displayed inequalities; adding w>=p-k1 gives exactly rL*p<=k1 and p<=k0+dU. For p>0 the reconstruction r=1-w/p,delta=w-D,b=k0-p attains every retained candidate. At p=0, feasible intervals force w=0 and delta=-D; choose any permitted r. This covers zero coefficients, empty sets and singleton boundaries without dividing by zero.

Sharpness means a witness exists for each retained p; it does not mean one common witness attains all p simultaneously. Correlated nuisance restrictions or a residual activity shared across genes require their joint constraints. Current text states the independent-range assumption, separates physical nonnegativity and warns against marginal-to-joint inference. No empirical nuisance bounds or confidence intervals are asserted. The reciprocal contrast is invariant only under one common positive lifetime-scale multiplier, and source fit/growth/uncertainty calibration remains unavailable.

## Strongest defensible claim

For these fixed algorithms and this operational reference-eligible population, residualized measured association gives a very small reduction in MAE predicting deposited depletion lifetime on1,112 model-heldout genes, and fails the declared practical criterion. It does not establish a calibrated causal contribution, true-isoform adjustment, a folding mechanism, biological-replicate inference or experimental replication. The failure also does not refute CNOT3 biology or prove the absence of nonparametric conditional information. The present manuscript, supplement and complete report respect these boundaries.

## Limits and closure

Only existing authorized phase outputs, manuscript/report sources, user-facing CSVs and prior source-audit summaries were read. The CSV audit includes already reported gene-level outputs and the supplied49-row arithmetic table. No original XLSX, new source measurement, model, outcome selection, fit or numerical hypothesis search was opened or run. Universal kinetic conclusions were checked by direct interval elimination; existing finite test receipts were not substituted for proof.

Review hashes for the current corrected TeX, abstract/discussion, supplement, both Markdown reports, JSON evidence and CSVs are in `checks.json`, with pre-correction manuscript hashes preserved. Final archive/version availability text is a release-stage item: the existing data-availability paragraph still names the older version1.0/e9846dc archive and should not imply that it already contains this later CNOT3 checkpoint. No archive/publication action was performed. Stop here; no additional biological expansion or tuning is warranted by this result review.
