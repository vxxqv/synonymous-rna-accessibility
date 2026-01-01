# FAS manuscript transfer and final figure review

Verdict: PASS for the mathematical transfer, reported counts, and corrected final figure. No material proof, arithmetic or visual defect was found. Two small interpretation clarifications are recommended below; neither changes the mathematics or numerical result.

Reviewed manuscript SHA256: `d5f2a59aabcf52fe74411dfb0c9ad3abd1bb4ab081debce35ad3d2aed4ad96dc`.

## Findings first

1. Lines159-162: “independent-coordinate envelope” correctly describes allowing each vertex coordinate to vary separately, but could be mistaken for statistical independence. Prefer “Cartesian sensitivity box, the product of the vertex intervals, without a statistical independence assumption.” The following explanation that it contains combinations not observed together is correct. The observations themselves remain dependent scenarios, and no confidence guarantee is asserted.

2. Lines191-192: “the strongest alternative” ranks a source-recovery explanation without having quantified its relative support against all alternatives. “A key unresolved alternative” would make the epistemic limit precise. This does not change the existing valid warning that the ordinal contradiction cannot identify splicing-specific mechanisms.

3. Optional proof completeness, after line110: the finite nondecreasing response assignment can be extended to a function on all real latent scores by constant tails and stepwise or linear monotone interpolation. The manuscript already constructs everything needed, so this is not a proof gap.

## Mathematical checks

- Lines87-113: the null permits an arbitrary **common nondecreasing** h, including plateaus. For each forced strict envelope separation U_i<L_j, latent order must be strict. Finitely many such margins can be rescaled to unit margin; this is not an effect-size calibration. Conversely, prefix maximum of finite lower endpoints over preceding/tied latent scores lies inside every interval: if it exceeded U_i, the corresponding forced order would reverse the assumed latent ordering. Latent ties give equal assigned responses. Infinite upper endpoints impose no outgoing forced order. No-forced-order feasibility includes a common constant response. Thus the necessity/sufficiency transfer is valid for the finite domain with valid closed intervals and finite lower bounds used here.
- Lines116-138: the five feature-row differences cancel exactly; the five required latent differences would all be positive. Equal weights1/5 give the normalized exact dual contradiction. All listed ratios, the 0.0075674 absolute gap, vertices42/51 and the rounded1.36% relative separation agree with the exact certificate. These margins remain observed separations, not uncertainty bounds or optimized robustness radii.
- Lines140-165: lower and upper field capacities bound every field in the Cartesian box. A guaranteed route exists when tau<=C_L. A guaranteed restriction with endpoints accepted throughout the box requires C_U<tau<=min(L_s,L_t). Weak threshold acceptance and the strict open lower endpoint of the restriction interval are correct, including equality. Outside those conditions the text retains uncertainty or lack of guaranteed endpoint acceptance. Zero primary certified restrictions is not universal biological connectivity. Secondary scenario barriers and a point barrier do not establish that a particular obstruction requires nonadditivity.
- The scope is the full72 native-reference translation class within the3072 design, not the full synonymous graph. The two P3 gates allow scalar unary valleys, so the wording correctly avoids attributing full-domain disconnection to intercodon interactions. The nucleotide/codon additive span equivalence is specific to the three-corner factor geometry and is correctly identified as a check rather than an additional biological model.

## Counts and source-scope checks

The transferred values agree with the saved independent exact audit:72 vertices,204 edges,2556 endpoint pairs,74 primary forced orders, five supported orders, zero primary envelope-certified barriers,27 dependent scenarios with62 to307 barriers,73 point barriers,1944 ratios,648 per-output relative scores,72 means,144 finite endpoints,30 dual certificates and74124 pair records. The saved source receipts have3071 rows in BR4,BR5,BR7,BR9,3072 in the other files; all72 saved selected variants have positive counts in all twelve files. Saved whole-file and TR1-aligned totals agree. This checks manuscript transfer from saved evidence, not a new raw count-file extraction or historical author-workflow authentication.

The text retains the headerless-import qualification, adaptive retrospective source selection, shared preparations, empirical ranges rather than confidence intervals, lack of independent biological validation, no calibrated PSI/protein function and unmeasured bypasses. Claims about closest prior studies and the original downstream-import behavior were not independently re-retrieved in this bounded mathematical pass; they rest on the separate source audits.

## Final figure delta

Reviewed PDF SHA256 `07f7248e83b46bea0b399133df23fcc37e64cd1f24f3b50f55b8f74714dfcc8c`; PNG SHA256 `686e82a98bf9db9316c2c703411c62f09b0d1f05e1437989915c033de5edccd9`; renderer SHA256 `12cb7095676f4d041038dfd8be7b90c494a3044cea3c5604304d0f439d40a74e`.

The full PDF was rendered to2055 pixels wide, corresponding to174mm at300dpi, and visually inspected. The corrected scope title/subtitle, feature-row wording and shortened lower-right line are present. No note/title overlap, clipped glyph or cropped extremum remains. All text boxes lie inside the612 by662.4 point media box: x2.084 to601.825, y12.488 to653.926. Smallest native text is9 points, or7.253 points at174mm; this matches the declared final typography choice, although it is smaller than an8-point preference. No8-point guarantee is claimed.

All four numerical CSVs remain byte-identical to those in the initial figure audit. The provenance source hash matches the authorized result and the provenance PDF hash matches this final PDF. All72 variants and all27 heatmap scenarios remain included. The five comparisons are the supplied exact certificate, not an optimized/minimal-effect selection. Caption labels match panel symbols and correctly define percentages as100(L_j/U_i-1), not uncertainty. The earlier clipping finding and receipts remain unchanged in the figure-review history.

This pass inspects the standalone figure and TeX source. It does not substitute for visual QA of the subsequently compiled manuscript page or formal revalidation of unrelated theorems.
