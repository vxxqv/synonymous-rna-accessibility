# Actual COMT saved-RF trace audit

**PASS for saved-state prediction and the declared source masks.** The sealed reviewer code was run once after root's explicit actual RF-result authorization. No forest was trained or tuned, no source was reacquired, and no predictive assay model or normalization repair was run.

| Check | ROI 1 | ROI 2 |
| --- | ---: | ---: |
| Extracted rows checked | 198,169 | 129,415 |
| Eligible prediction rows checked | 198,169 | 129,415 |
| Maximum absolute probability difference | 0 | 0 |
| Class mismatches | 0 | 0 |
| Rejected SNP identities and counts matched | 6 | 1 |
| Source frequency-pass rows matched | 155,732 | 97,365 |
| Saved training-row metadata | 629 | 783 |
| Saved seed metadata | 9 | 9 |
| Trees | 500 | 500 |
| Saved/aggregated mtry | 3 | 3 |

All 327,584 extracted record identities and order, source raw variant types, independently encoded feature columns, saved probabilities, >=0.5 classes, full frequency masks and global rejection masks agree. The independently rebuilt gDNA false-call aggregation agrees with both rejection files: only per-sample singleton SNP false calls contribute, with >=4 among the seven declared ROI 1 gDNA libraries and >=2 among the four ROI 2 libraries. These row counts include dependent per-base/allele/sample records; they are not independent biological observations or synonymous-target counts.

The probability tolerance was frozen at absolute 1e-14 before access; no difference required that tolerance in this run. Class labels, masks, identity sets/counts and feature columns were required to match exactly. The reviewer constructed factors/features without importing the author preparation function or RF adapter, then used the same pinned randomForest engine on the saved forest objects. This does not independently implement the forest prediction algorithm or validate training/tuning, classifier accuracy, historical equivalence or biological calibration. Seed9 is checked as saved metadata, not authenticated RNG history.

## Provenance and exposure

The extraction seal matches root's `4e365b305c356a25762c540092083c1a6e87bbc89041c97eb90a2e8f0436d7f0`; the RF seal matches `091b66e102669f466121244e65a85c5434eafd514e903285074515a51bdd7564`. The common protocol, amendment, integration addendum, contract and code pins agree across phases. RF points to that exact extraction, and both stages point to the frozen runtime. The saved models, every sealed RF/extraction payload and runtime package/executable hashes were verified before replay. `actual_authorization_v1.json` records scoped root authority and pre-access verifier identity. Public revision `13b900f704b76b1b3d3059a3bbec8ac1f510b715` is root-supplied and was not independently reauthenticated by Git.

The unchanged pre-access seal is `4e369a6b7a61d517e6af1d4e544bf05ed8ef9980572401fa77e8024eba506dfb`. Python verifier SHA256 is `04540a406360b55bbec3d1b80486020f98904520015a353f436a919f1911eff4`; R verifier is `f3cefa6185e274030f69849f8c629cfb0897c839d91c3a2777ea14a5f17a52c2`. No method, tolerance, code or case selection changed after actual access. `actual_v1/verification.json` is `669d7464ed71ff96be8ba1f2406ff78f09fd4367314c722ba2594f96288075df`.

The authorized extraction tables and saved RF objects/outputs were read. Raw archive bodies, original source simulations/truth files, source score tables, outcome missingness, later contrasts, prediction/evaluation artifacts and the normalization failure's row values were not opened for this audit. Classifier CAF/covariates and saved probabilities/masks are now exposed within this authorized scope. Reading the extracted tables as text does not make the rest of this source a held-out biological validation resource.

## Warnings and limitations

The fresh R process emitted four startup warnings that `LC_COLLATE`, `LC_CTYPE`, `LC_MONETARY` and `LC_TIME` could not be set to inherited `C.UTF-8`. They are retained verbatim in `actual_v1/execution.log`. The verifier's subsequent explicit `LC_COLLATE=C` call succeeded, and all feature columns and prediction outputs matched exactly. No warning was hidden, retried away or used to revise a setting. This trace does not establish a fully identical OS locale/base-R environment.

The later unique post-collapse-DP normalization assumption failed in production, according to root's message. This report does not inspect, repair or reinterpret that failure, and RF trace agreement cannot turn it into a completed assay-prediction analysis. Root must preserve and separately review that gate. No additional actual run is needed for this RF verification.

Inspect `actual_v1/trace/trace_1.tsv` and `trace_2.tsv` for every retained probability/class comparison, `summary.tsv` for full counts, `verification.json` for input/output hashes and `after_check.json` for post-run immutability checks. The reviewer stops here, with no commits, publication or new biological comparisons.
