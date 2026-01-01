# Figure 12 handoff

Status: ready for root visual review and integration. No blocking issue remains. No manuscript, repository, model or outcome source was edited. All work is under `work/decay_cycle4/figure12/`.

## Authoritative files

- Portable renderer: `render_figure12.py`.
- Frozen panel data: `panel_data/Figure_12_v5.json`.
- Final exports: `final_v5/Figure_12.pdf`, `final_v5/Figure_12.svg`, `final_v5/Figure_12.png`.
- Rendering receipt: `final_v5/render_receipt.json`.
- Final automated layout/replay receipt: `qa_v5/layout_receipt.json`.
- Actual PDF preview: `qa_v5/Figure_12_pdf_review.png`.
- Grayscale preview: `qa_v5/Figure_12_grayscale.png`.
- Full source and export hashes: `integration_manifest.json`.

The earlier `exports_v2`, `exports_v3`, `final`, `qa` and candidate panel JSON files are layout history only. Do not integrate them. The original and final panel dictionaries compare exactly, so layout changes did not change any scientific value.

## Caption draft

**Figure 12. A small exploratory prediction gain does not identify a pathway contribution.** (a) Equal-gene mean absolute error (MAE) of natural-log CNOT3-knockout half-life for all nine frozen models on 1,112 model-held-out genes in 1,109 protein-similarity groups. Lower MAE is better. Filled gray markers denote the five baseline controls, teal denotes the primary association model, and open markers denote the three frozen sensitivity variants. The primary model fails the practical criterion requiring at least a 5% relative MAE reduction over each required baseline. Its reduction is 0.113% against the codon/core baseline and 0.042% against the amino-acid/core baseline. (b) Relative MAE reductions against the matched codon baseline are 0.113%, 0.100% and 0.319% for the primary, no-SE and nonlinear-IP association variants, respectively. The dashed 5% line is the primary practical criterion shown for reference, not a newly assigned success criterion for the sensitivities. (c) Rank-percentile comparison of lifetime change, Delta = log2(t1/t0), and the signed dimensionless reciprocal contrast, R = ln(2) h_ref (1/t0 - 1/t1), for all 5,416 eligible genes, with no point excluded. Here h_ref = 5.17814182006782 in the same time units as the half-lives; exact ties receive average ranks, divided by 5,416. Spearman rho = 0.826797. These quantities are algebraically coupled through R = ln(2) h_ref/t0 [1 - 2^(-Delta)]; the rank plot does not display their magnitudes or establish a mechanism. (d) A reviewed synthetic counterexample under k0 = b + p and k1 = b + delta + r p. Even with r = 1/2 fixed and nonnegative background rates, the identical observed total rates (k0,k1) = (3,2) allow p = 2 or p = 3 when compensation delta is not fixed. Rates use arbitrary consistent units. This is a mathematical example, not a fitted CNOT3 decomposition. The empirical comparison is exploratory and source-exposed: model-held-out genes do not imply historically untouched outcomes. No calibrated biological uncertainty, confidence intervals or causal pathway attribution is supplied. The failed original source-reproduction gate is retained separately.

## Numerical and visual checks

- Frozen source hashes were checked before and after extraction. Only the authorized saved evaluation, independent verification, cohort receipt and previously reviewed kinetic derivation were used. No assay workbook, target-table source, new fold or model fit was opened or run.
- Independent verification status is passed and covers this exact evaluation and prediction hash. The figure does not independently authenticate historical Git access order.
- All nine MAEs and all three matched relative gains come from the frozen evaluation. Gain arithmetic was checked independently during panel extraction. Primary success remains false.
- The scatter includes every one of 5,416 finite, unique eligible GeneIDs. Average-rank correlation reproduces the frozen rho within 1e-12. The coupling identity was checked on every saved row. Raw Delta range is [-2.8864013573136056, 4.496929869997562]; raw R range is [-0.8091713985061316, 1.7766805636246232]. Plot bounds [-2,102] on both percentile axes omit no point or outlier.
- Both synthetic decompositions were verified using rational arithmetic, including b >= 0, p >= 0 and b + delta >= 0.
- The final PDF is one page, 493.228 by 464.882 points, equivalent to 174 by 164 mm. PNG dimensions are 4,178 by 3,938 pixels at 610 dpi. The PDF and SVG contain vector plot elements; the 5,416 scatter markers are not rasterized.
- All 73 visible text items are at least 8.5 pt. Automated text bounds are inside the page, and no pair overlaps by more than 0.6 pixel in both dimensions. These are text-geometry checks, not a proof of universal renderer behavior.
- The actual PDF was rendered by Poppler at 180 dpi and inspected in full, as was the grayscale image. Panel letters, equations, marker/value columns, tick labels, threshold line, synthetic table and caption-adjacent footer are legible and separated. No remaining clipping, line/text overlap or color-only distinction was observed. No claim of theorem revalidation is made by visual inspection.
- A fresh render produced byte-identical PDF, SVG and PNG assets. A failed earlier geometry check was corrected before finalization and remains in `qa/` as historical evidence.
- `pdffonts` is unavailable locally. The renderer explicitly uses PDF font type 42 and SVG text, but no independent font-embedding utility result is claimed.

## Portable render command

Python 3.12.14, NumPy 2.5.4 and Matplotlib 3.11.2 were used. The renderer needs only preinstalled NumPy/Matplotlib and the panel JSON. It performs no network or Git operation and does not use the original workspace layout. Use a fresh output directory:

```text
python render_figure12.py --data Figure_12_v5.json --expected-data-sha256 3f40c1d627effd1fb0934b6ff676952bcca7030ae88e00b774b5576d061dbf31 --out NEW_DIRECTORY
```

The renderer refuses changed data, changed renderer identity, changed frozen source pins, an existing output directory and disabled assertions. The supplied JSON and renderer must be copied without changing their bytes. The preparer and QA helper are provenance utilities, separate from the portable rendering entrypoint.
