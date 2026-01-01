# D1 draft: nested local-descriptor order models

Status: design-only proposal. Root review, implementation review and an explicit new freeze are required before any fit. These source outcomes were already exposed in B1. A later run is exploratory reuse, never prospective validation. No source values have been read in this design audit.

## Fixed units and inputs

Retain all 66 selected local states, all 30 intended-edit groups, exons 17, 25, 27, 39 and 42, both K562 and KBM7, and all original B1 interval/missingness rules. No new selection or score recomputation. Read the independently verified saved B1 exact interval records through an explicit input path and SHA-256 pin to be supplied by root. The protocol is not executable until that result path, hash and result schema are authenticated. Do not reopen the phenotype workbook.

The matrix rows and numerical entries are frozen in `inventory_v1/nested_matrices.json`. Their origin is `exact_v1/mapping.json` and `exact_v1/certificates.json`; file hashes are in the receipt and final manifest. Row order is K562 then KBM7, each with lexicographic canonical state IDs. Every local base outside its counted DNA word remains unknown. The model describes the two recorded allele descriptors; it does not assert complete haplotypes, invariant protein, or absence of other changes.

Use the four existing widening factors exactly: gamma in {1, 5/4, 3/2, 2}. A finite positive interval [L,U] becomes [L/gamma,U*gamma]. Unknown remains (0,infinity). Every interval bound and decision uses exact rational arithmetic. Retain all states and failures, including unknown states. No additional tolerance, selected threshold or outcome-dependent row removal.

## Actual nested model classes

For each exon, let X have an intercept and the recorded physical allele indicators. Exon-specific parameters are unrestricted. Separate nondecreasing response functions h_(cell,exon) map latent scores into that context's response intervals.

M0, shared descriptor additive: z_(cell,v)=X_v beta_exon, with the same beta in both cells. M1, context descriptor additive: z_(cell,v)=X_v beta_(cell,exon). M2, context-specific observed-pair model: augment X by the product of the two descriptor indicators for every observed physical allele pair, separately in each cell. Include every pair column, not an outcome-selected term. Each pair column is the indicator of one observed row, so M2 is a saturated latent-score reference, not an identified interaction coefficient model.

M0 is contained in M1 by beta_K562=beta_KBM7. M1 is contained in M2 by setting all added pair coefficients to zero. M0's column ranks are 25,4,10,9,17 by exon; M1 doubles each rank; M2 has twice the exon state count. Thus M1 and M2 have identical fitted latent spaces for exons 25,27,39,42. Only exon17 can distinguish them, through one constraint per cell. A full-rank single-context design can represent every finite local order, so it cannot reject the interval-order model when its intervals are nonempty. This expected non-rejection is a design limit, not evidence of additivity.

Cell-specific intercepts are absorbed by the independently unrestricted h functions and are not a distinct ordinal hypothesis. Likewise, adding a common constant within an exon has no effect on its strict order constraints. No cross-exon order is tested.

The previously computed unrestricted common-order comparison remains a separate, weaker reference: M0 implies a common order, but the converse is false. Use all five exons and all four gamma values in any comparison; do not replace M0 with that prior test.

## Exact equivalence to strict linear feasibility

For finitely many nonempty intervals I_i=[L_i,U_i] or (0,infinity), define a forced pair i>j precisely when L_i>U_j, with infinity never a finite upper bound. If responses y_i=h(z_i) lie in the intervals and h is nondecreasing, then a forced pair requires z_i>z_j. Conversely, suppose every such forced pair has z_i>z_j. Order the distinct latent levels, group equal scores, and choose a common positive response at each level by taking the maximum of its lower bounds and the previous chosen value. Start from a positive rational delta smaller than all finite positive bounds, or delta=1 if none exist. A violated upper bound would imply a forced pair whose earlier or tied score is not larger, a contradiction. Ties within a level therefore have a nonempty common interval; unknown open-zero intervals receive a strictly positive value. Extending these finitely many levels by constant pieces gives a nondecreasing h. This establishes equivalence, including ties and unknown intervals.

For each forced pair, form the exact row a_ij = design_i - design_j and require a_ij beta > 0. There are finitely many pairs. If any strict feasible beta exists, the minimum positive gap is positive and scaling beta by its reciprocal yields a_ij beta >= 1. Conversely the unit-margin system implies every strict inequality. No bound or regularizer on beta is imposed. The unit margin is a harmless mathematical normalization, not a biological effect size. With no forced pairs, beta=0 and the monotone completion above are a valid witness.

Rational matrix entries imply that a nonempty unit-margin polyhedron has a rational point. Every feasible classification must save beta, all exact latent values, every forced-pair gap and the completed responses in each context. Every infeasible classification must save a nonnegative rational vector lambda with sum(lambda)=1 and sum(lambda*a)=0. Then multiplying A beta >= 1 by lambda yields 0>=1, an exact contradiction. A floating solver status alone is insufficient. Certificate reconstruction failure is an unresolved computational result, not infeasibility.

## Complete reporting and failure rules

Report every exon x model x gamma record: 5 x 3 x 4 = 60 joint-two-context comparisons, plus the two individual-context feasibility records for each M1/M2 comparison. Do not sum these as independent experiments. The global all-exon result is the conjunction of the per-exon results because all parameter blocks and response functions are exon-specific. Its certificate can use one failing exon or concatenate all feasible witnesses.

M2 should always be feasible for nonempty intervals by its row-indicator columns. Any contrary result is an implementation or input error that stops interpretation. M1 should always be feasible in the four saturated exons. A shared-model rejection can reflect context-specific ordering or calibration/measurement dependence. An M1 rejection in exon17 excludes only the specified local-descriptor additive latent model at that widening level. It does not identify a physical nucleotide interaction, folding mechanism, protein consequence or independent biological cause. M2 feasibility alone establishes nothing about interactions.

Do not interpret a raw observed-scale balanced contrast as a latent interaction: a strictly monotone nonlinear response can produce it from an additive latent field. Raw contrast calculations, if separately proposed later, require their own explicitly descriptive endpoint. This protocol tests only compatibility with interval orders.

Mandatory checks before any run: exact synthetic witnesses/duals; row and column permutations; all missing/tied fields; expected saturation; reference coding span invariance; strict-to-unit-margin normalization; gamma nesting; all 66 metadata identities; saved interval hash; no overwrite; and independent certificate verification. Increasing gamma only removes forced pairs, so an earlier feasible model cannot become infeasible later. Preserve every failure without changing gamma, model or sample set.

The finite algebra is standard order feasibility, linear-model estimability and Farkas certification. The proposed increment is the complete, explicitly limited application to this authenticated selected design. It is not a new general theorem.
