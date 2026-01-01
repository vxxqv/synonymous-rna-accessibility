# FAS figure and CSV audit

Verdict: numerical PASS and current-layout PASS; two scope-label clarifications and final-size typography adjustment recommended before release. No source files were edited. These findings refer to the hashes below and not to a subsequent regeneration.

Inspected the original PNG and a fresh 2055-pixel-wide PDF render corresponding to 174 mm at 300 dpi. No clipping, text/marker overlap or missing extrema was observed. PDF text boxes are inside its 612 by 662.4 point media box (x from 2.084 to 605.580, y from 12.488 to 653.760 points). Panel a includes all 72 variants, all exact range endpoints and the correct native reference at vertex 46 and score 1. No outlier was cropped. Teal points, gray ranges, navy squares and red gaps remain distinguishable; heatmap cells all have numeric labels. The missing colorbar is harmless because counts are explicitly printed and the common color scale is fixed from 0 to 320.

## Actionable wording and size findings

1. Panel b title, currently “Five orders exclude the monotone unary model”: explicitly name the **common monotone codon-additive model**. For example, “Five orders exclude common monotone codon additivity.” The exact conclusion requires a common response h applied to a sum of codon effects; arbitrary vertex-specific unary scores are not rejected. This is a scope clarification, not a numerical defect.

2. Subtitle, currently “Complete native included-protein class”: qualify this as **within the published design**. These are all 72 native-protein-preserving sequences in the authenticated 3072-sequence design, not every possible synonymous sequence in the ambient coding space. Suggested wording: “All native-protein variants in the published design | minigene RNA-enrichment assay.” Actual reporter protein expression is not established.

3. At the requested 174 mm final width, the 8.5 inch canvas is scaled by 0.805929. The 8.3 point footer becomes 6.689 points; 8.5 point notes/legend become 6.850 points. They are legible when magnified, but small at final print size. Native size 9.93 points is needed for a final 8 points. Increase the small native fonts or author at final width and inspect the resulting layout again. Existing margins have no collision; increasing font size may require line wrapping.

4. Optional precision: “coefficient difference” in panel b note is more clearly “feature-row difference.” The cancellation is of the five X_j-X_i rows; model coefficients b are unknown. Their equal weights make both weighted and unweighted row sums zero, so the current sum statement is mathematically valid.

## Exact numerical checks

All four CSVs agree with the independently verified result: 72 measurement rows, all 2556 primary endpoint pairs, all 28 scenario/point summary rows, and all five certificate orders. Exact fractions, saved count vectors, all scenario ratios, source point means, dual differences and weights agree. Complete row ordering is preserved. The figure's source hash equals the authorized result hash, and its provenance PDF hash matches the inspected PDF.

The 1.36% label is exactly the rounded relative excess

`100 * (L51/U42 - 1) = 1.3645978866369182%`.

The separation ratio is `1654983897462/1632704052467 = 1.0136459788663692`; the absolute endpoint-envelope gap is about 0.007567413977870934. None is a calibrated uncertainty bound, optimized robustness radius, interaction magnitude or confidence statement. Panel b's note correctly distinguishes separation ratios from confidence. All five percentage labels use the same relative-excess convention.

The primary envelope still has zero certified barrier pairs. The heatmap includes every one of the 27 dependent scenarios; its values range from 62 to 307. The point-summary count is 73. The figure does not average these into an independent replication claim or claim that null rejection establishes a barrier. The footer correctly retains the lack of calibrated confidence intervals and independent biological validation.

## Exact inspected hashes

- PDF: `0f228cfe4581e80448fc8602a4c7abb90452ea07741675d84b1127b31d2d8597`
- PNG: `9cf66cb599a18918f8a2c1b951add5a67e13a95421811b91f4eee159e72368f8`
- Renderer: `77444ad59b92eab28aabd6baa8ce9e40f9f602b6d05ac33d7e313c775e4cf1e8`
- Original result: `3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466`

Full CSV/provenance hashes and size calculations are in `checks.json`. This audit checks visual rendering and faithful presentation of the already verified finite result. It does not revalidate molecular mechanisms, source extraction, historical access order or biological uncertainty.
