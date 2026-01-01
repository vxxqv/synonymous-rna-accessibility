# Independent source discrepancy audit

The frozen failure is reproduced. The independent reader finds the same5,480 cohort GeneIDs and the same49 stored-ratio tolerance failures. No implementation error explaining either failure was found. A10-significant-digit intermediate calculation accounts numerically for every stored ratio, but the published5,481-gene count remains unauthenticated. No scientific rule or failed gate was changed.

## Count:5,480 is not an implementation off-by-one

Within the frozen exposure ledger's7,359 finite-IP candidate rows,5,480 have positive finite WT and KO half-lives and a finite stored log2 ratio. The other1,879 have **all three B/C/D cells absent**. There is no single candidate with only one half-life or only the stored ratio missing. All7,359 candidates independently satisfy finite I/J/K, and their archived STAT agrees with IP log2 enrichment divided by SE under the original tolerance. The independently extracted cohort and its raw numeric export match the frozen reproduction exactly.

The exposure ledger's paired selector contains5,481 row numbers because it deliberately includes header row1 and5,480 data rows. The reproduction correctly excludes that header. This explains that ledger count; it is not evidence that the publication included a header. The article's Figure2 caption reports1,000 high-ranked genes plus4,481 remaining genes. No published5,481-ID membership list or exact original missingness/read-filter implementation is available in the inspected evidence, so the identity of a supposed extra gene cannot be established. A caption typo, a different original universe and an unrecorded source step remain distinguishable possibilities without evidence selecting one. No alternative cohort was searched or constructed.

This audit independently rechecks the already exposed finite-IP candidate universe. It does not semantically decode I/J/K outside that frozen selector or claim a new independent search across the remaining workbook rows.

## Ratio arithmetic and rounding evidence

The workbook documentation defines D as log2(KO half-life / WT half-life). An independent SAX parser retained original XML numeric strings.80-digit Decimal arithmetic reproduces all49 failed IDs under the unchanged tolerance `1e-10 + 1e-8*max(abs(computed),abs(stored))`. The frozen author's computed ratios differ from the high-precision calculation by at most8.92e-16, so binary logarithm arithmetic is not an explanation for the larger source residuals.

- Maximum absolute stored-minus-recomputed residual across all5,480 rows:1.229557808497134e-9.
- Among the49 failures, absolute residuals range from1.320991054495693e-10 to9.507555238026865e-10; the maximum tolerance exceedance is4.43549-fold.
- All49 have absolute stored log2 change at most0.05963903. Their near-zero values explain why the frozen relative-plus-absolute tolerance detects these small residuals.

All paired B/C/D cells use Excel number format **General**. This supplies no fixed display-decimal or10-significant-digit precision rule. Conditional half-unit bounds at the *last written raw XML decimal* cannot reconcile any of the49 failures. Therefore ordinary rounding at the retained15-or-so-digit XML precision is insufficient; a more substantial intermediate quantization would be needed. These bounds test an explicitly stated rounding hypothesis, not measurement uncertainty.

Root subsequently supplied a specific post-exposure diagnostic: round each raw half-life to10 significant digits, then calculate log2 of their ratio. Independently applying decimal `ROUND_HALF_EVEN` reproduces the stored ratios for **all5,480 rows**, with maximum absolute difference5.286335695964195e-15. The diagnostic agreement threshold was2e-13. Every failed row's raw values,10-digit values, high-precision recomputations and residuals are retained in `observed_v2/ratio_mismatches.csv`.

This is strong numerical compatibility with a10-significant-digit intermediate. It does not prove the unpublished generating algorithm, identify the software/export stage, authenticate unrounded half-life estimates or make the original frozen gate pass. The article says its algorithms are available on request; this audit did not obtain them or contact anyone. The rounding hypothesis was examined after exposure and cannot retrospectively become the original prespecified arithmetic rule.

## Scope, artifacts and stop

Input workbook SHA256: `a338e7211fcb0ad6bed0b1bfdb57433a5a43266f3a8d93e473b1bf50c583d712`; scientific freeze commit: `cc9cb18954026b9fe6199afe9820920823aaa767`. Only A–D and I–K on already exposed candidate rows were decoded. F/G were unnecessary and remained unopened by this audit; E,H,L–X and all other endpoint rows remained undecoded. Worksheet/style metadata and only needed shared strings were read. No author parsing or analysis functions were imported.

`audit_source.py --out NEW_OWNED_DIRECTORY` reproduces the bounded check. The latest receipt is `observed_v2/summary.json`; the49-row CSV and candidate-exclusion ledger are beside it. The initial script/output remain historical. Their numerical CSV and exclusion ledger are byte-identical to v2; v2 adds the explicit output argument and native-versus-Decimal comparison. Five synthetic parser/arithmetic tests supplement source inspection. Inputs were unchanged before/after both runs.

No modeling, grouping, prediction, threshold adjustment, source replacement or scientific protocol edit was performed. The appropriate status remains **held for source diagnosis**. The predictive hypothesis has not been tested, and this source-gate failure is not evidence against that biological hypothesis. The unresolved count prevents closing the prespecified reproduction gate even if root separately records the rounding-compatible ratio explanation.
