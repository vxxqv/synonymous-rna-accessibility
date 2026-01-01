# Independent completed-cohort review

**Verdict: passed for the declared operational cohort and exact grouping rule.** The independent calculation completed in24.41seconds using only Python standard-library logic and importing none of the author's cohort or prediction functions. No material discrepancy was found. This review used the final completed grouping, not its last in-progress checkpoint.

| Check | Independently recovered result |
| --- | --- |
| Operational source and original-ID intersection | P=5480, S=5480 |
| Reference eligible cohort | A=5416, A/S=0.9883211678832117 |
| Exclusions |38 selected representatives not protein_coding;14 incomplete UTR annotations;12 missing selected references |
| Training / holdout genes |4304 /1112 |
| Training / holdout groups |4280 /1109 |
| Inner assessment genes, folds0 through4 |855,869,875,882,823 |
| Predictor schema and content |All92 fields of every row agree with the allowed source fields and pinned reference; no KO/change-score field; zero-length UTR GC remains missing |
| Protein/reference checks |5415 unique proteins;27080 transcript/CDS/UTR/protein hash comparisons |
| Candidate/edge verification |6502921 shared5-mer posting pairs; all72 candidates independently regenerated; all18083136 edit-distance DP cells independently calculated;33 qualifying near edges and1 exact-identity group |
| Connected components |5389, with every group ID and member set reproduced |

All source/retained/excluded identifiers and exclusion reasons were checked, including sequential counts. Both hash partitions were recalculated from all32 SHA256 bytes with the original salts. No gene or exact stored component crosses the outer split or its training fold. All gene/group/fold coverage guards passed independently. A check of the complete candidate set, rather than only accepted edges, supports the absence of a missed edge under the specified full-length edit-distance threshold.

The independent distance routine uses the elementary full Wagner-Fischer recurrence, not the native aligner used by the original implementation. Candidate construction uses independent dense shared-word counts with the established lossless5-mer lower bound. For distinct strings capable of passing the10% rule, maximum length is at least10 and `maximum_length-4-5*floor(maximum_length/10)>0`, so a passing pair must share a5-mer. Exact identities are handled separately. Connected components were rebuilt by graph traversal, separately from the author's union-find implementation.

Seven actual scientific input files and six reviewed final output files were byte-hashed against the completed immutable receipt. The source gate and independent source-audit semantics were handled in the separately completed gate/source reviews, rather than reread here. Targets, model predictions and evaluation outputs were not opened.

The flat source-row reader explicitly skips the HL_CNOT3KO literal without numerical decoding. Its poisoned-literal fixture confirms that this field does not enter its numerical parser. All5480 KO literals were skipped; zero KO values were decoded and zero target files were opened. Only allowed baseline/enrichment fields were compared with predictor files. Positive paired-half-life validity remains supported by the separate source audit, not a new outcome read in this review. No model was fitted.

This pass does not resolve the missing historical5481st gene, source half-life uncertainty/calibration, broader paralog/domain relationships, expressed isoforms or experimental dependence. It verifies the stated operational cohort and computational grouping, not biological replicate independence or independent experimental validation. The reviewer had inherited project context and authored only the operational guard adaptation; the underlying eligibility/grouping/reference implementation came from another lane.

Executable evidence: `check_cohort.py`; exact candidate calculations: `independent_candidate_distances.json`; machine-readable counts, source/output hashes and exposure record: `result.json`. These are the only new calculations in this review. No original code, source input, manuscript, repository or publication was changed.
