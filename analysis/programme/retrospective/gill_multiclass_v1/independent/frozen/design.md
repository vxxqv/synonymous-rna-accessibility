# Independent Gill result-verification contract

This contract is fixed before this reviewer opens RNA D, exported validation counts, author class metrics or the author validation summary. Root has reported that its validation execution completed, but has not supplied numerical results to this reviewer. The reviewer has inspected the committed protocol/inventory and schemas/receipt of the permitted training exports. No author scientific implementation is imported. This is separately authored with inherited project context, not a blinded external review.

## Pinned inputs and membership

Protocol SHA256: `4ab4e85c14203c0ab8c983e8433c766f1566f373a335fd3fda4fb3981371c1aa`. Decompressed inventory SHA256: `2e00b698dee109e357dc01435d5eacf8b341388f36463970454ade230664f350`. Training receipt SHA256: `c6fb3309a93635efab38f1469516617c3f1b4f857c88b4fbfff91be21115bdef`. Source workbook identity is inherited from the frozen protocol/receipt; this verifier will not open that workbook.

The full fixed membership is2176 classes:1124 development/20028vertices and1052 validation/18310vertices. All38338 keys must occur exactly once in the selected class inventory, with complete declared class membership and exact training/validation exported key sets. The training export must have B,C,E,F,G only; validation must have D,G only and reproduce the training export's G. All selected counts must be nonnegative integers; zeros remain. Combined DNA E+F+G must exceed50. Missing/invalid rows or changed hashes stop verification rather than remove classes.

All selected keys are checked as complete Cartesian products of the four codon-variable blocks `key[0]`, `key[1:4]`, `key[4:7]`, `key[7]`. Within this fixed construct, these identify the four varying codons because omitted flanking bases are fixed. This checks the numerical factorization, not a new independent sequencing/translation audit of physical molecules. Each selected graph is rebuilt from one-nucleotide Hamming edges and checked for connectivity and exact frozen edge count. Expected factor ranks and residual dimensions must agree; the reference-protein stratum VP remains quarantined.

## Independent projection and endpoint

For a complete four-factor product, construct an indicator matrix with an intercept and one omitted level per categorical codon factor. Its rank is `p=1+sum_j(s_j-1)`, and `q=n-p` must match inventory and be at least3. Fit centered numeric fields with `numpy.linalg.lstsq` and reconstruct P and Q from this indicator design. The author implementation is not used. The protocol's conditional-mean expression follows because the four centered main-effect subspaces are mutually orthogonal under the complete uniform product measure. Metadata/synthetic preparation compares these two representations on every selected class.

For c in{0.25,0.5,1}, construct all six log2(count+c) fields and use

`T=(logB+logC-logE-logF)/2`, `V=logD-logG`.

Apply the linear projector consistently; Q annihilates the intercept. The primary numerator is `sum_g dot(QTg,QVg)/n_g`, with analogous pooled squared norms. Divide only if both pooled energies are nonzero. This is the specified inverse-class-size weighting of vertex products, so every class contributes its mean product. It is not an independent-gene sampling calculation. Do not discard any class because its individual residual energy vanishes.

Compute class residual cosines, residual prediction MSE `mean((QT-QV)^2)`, zero-residual MSE `mean(QV^2)`, centered whole-field cosine and fixed-T MSE, constant-predictor MSE, and centered codon-main-effect prediction loss. Preserve negative prediction gains and undefined zero-denominator cases. No fitted shrinkage or threshold selection is added. Recompute all3 RNA-only and DNA-only residual pair cosines and each class's raw-RNA midrank agreement, explicitly retaining undefined ranks.

## Normalization identity, all scenarios and influence

Scenario order is B/E,B/F,B/G,C/E,C/F,C/G,D/E,D/F,D/G. For each class/pseudocount, compute all45 unordered scenario pairs including diagonals. Save the four signed raw inner-product terms, their direct total and numerical identity discrepancy. Also compute the averaged-training versus D/G decomposition. Aggregate each signed term with1/n_g weights. Cross terms are numerical contributions, not causal attribution fractions.

For each partition and pseudocount, group the primary numerator/energy contributions by local central-AA stratum. Recompute every leave-one-represented-stratum-out pooled cosine, its range, and the largest fraction of training and validation residual energy contributed by one stratum. The eligible selected data represent162 development and154 validation strata. The inventory's200/199 figures describe assignment of possible strata before eligibility, not these observed denominators. Influence calculations are dependent sensitivity analyses, not confidence intervals or independent biological resampling.

## Exhaustive graph verification

No architecture subsampling is needed. Selected classes have20 factor-size architectures, maximum72vertices, and sum n^3=33733816. For every2176class and every3pseudocount, compute capacities for all endpoint pairs by max-min Floyd matrix closure on the independently rebuilt nucleotide graph. Check all9 scenario fields and their codon-additive projections, plus T and PT:20 field/graph calculations per class/pseudocount,130560 in total. V and PV are already scenario8 in zero-based order and are also reported as held-out fields.

Initialize an edge's capacity to the minimum endpoint score, a diagonal to its vertex score, and a nonedge to negative infinity. For each intermediate vertex k apply `C_ij=max(C_ij,min(C_ik,C_kj))`. The maximum endpoint-relative barrier is `max_ij(min(F_i,F_j)-C_ij)`. Keep the argmax pair and per-class maxima. Compare clique-product and nonclique-factor classes separately without presuming nonclique codon-additive barriers vanish. All extrema can be independently derived from the complete saved class arrays. This is exhaustive finite graph checking, not a new universal theorem.

## Comparison policy fixed before author results

The scientific formulas and independent evaluator are frozen now. A later adapter may map author field names/order to these fixed quantities after root supplies the receipt and schema. It must not change membership, equations, pseudocounts, graph coverage or tolerances from observed agreement. The independent evaluation should be written and hashed before reading the author's numeric summary, then compared against all author class arrays and aggregates.

Memberships, key order mappings, integer ranks/dimensions and declared counts require exact agreement. Projection-vector comparison uses a norm budget `delta=256*eps*||raw_logs|| + 1e-10*||centered_logs||`. Dot-product/energy comparison uses propagated budgets of the form `delta*(||u||+||v||+2*delta)`. Barrier differences are compared in input score units with `1e-9*max(score_range,eps*max_abs_score,tiny)` per field. Ordinary well-conditioned cosine/ratio comparisons require agreement within1e-9, supplemented by a recorded propagated conditioning bound when a denominator is small. Near-numerical-null individual class cosines are reported as conditioning/undefined-status discrepancies; this is not permission to exclude them, change their contribution or alter the primary endpoint.

Every discrepancy must be recorded. Do not widen a tolerance after seeing results, erase an undefined value, round a failure across the R=0.5 target, or equate algebraic residual dimensions with biological error degrees of freedom. The primary operational pass/fail is the unrounded frozen comparison R>=0.5; a threshold-adjacent numerical ambiguity must be reported rather than promoted. Preserve the author's failed targets and all sensitivity directions if reproduced.

## Execution status and remaining limitations

Preparation is metadata/synthetic only. `preparation_frozen/manifest.json` is the authoritative successful final-source run. Earlier preparation artifacts are retained. In particular, `preparation_final` records a failed reviewer assertion that all assigned possible strata would be represented after eligibility. Its exact source snapshot is `verifier_failed_assumed_strata.py`. Metadata showed162/154 represented strata, so that assumption was corrected before verifier freeze and before validation access. This was not an outcome-driven cohort change. An earlier successful preparation run refers to an older reviewer source revision and is historical, not the final-source validation receipt.

The evaluate phase requires explicitly supplied validation receipt/count paths and expected hashes and checks the frozen verifier source hash. It is not authorized to discover a validation directory or read the workbook. It can verify the exported counts and downstream results, but cannot independently certify correct source-workbook extraction without a separate permitted source audit. Preparation checks passing do not verify any observed numerical biological result.

Regardless of agreement, conclusions remain descriptive repeatability of processed measurements from one truncated reporter and local protein strata. They do not establish independent biological replicate identity, calibrated noise/selection decomposition, causal molecular interactions, a biological noise ceiling, or independent-gene/study generalization. Log-scale residuals can arise from a nonlinear response to a latent additive field. No p-value or confidence interval will be introduced.
