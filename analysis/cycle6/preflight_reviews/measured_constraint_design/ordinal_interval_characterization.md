# Exact finite interval test for a monotone additive null

This is a self-contained application of finite order representation and linear alternatives, not a priority claim or new general measurement theorem. No outcomes were used.

## Contract

There are n>=1 rows X_i in Q^p, with unrestricted real coefficients b. Every interval is [L_i,U_i], where L_i is finite, U_i is finite or +infinity, and L_i<=U_i. The biological application has L_i>=0, but the proof only needs finite lower bounds. Seek one function h:R->R that is nondecreasing, possibly constant or with plateaus, such that h(X_i b) belongs to every interval. No coefficient box, regularization constraint or fixed sigmoid is imposed. Feature rows need not be distinct; p=0 is allowed.

Define E={(i,j): U_i is finite and U_i<L_j}. The exact criterion is feasibility of

`(X_j-X_i)b >= 1` for all (i,j) in E.

## Proof of necessity

If h(X_i b) is in its interval and U_i<L_j, then h(X_i b)<h(X_j b). Nondecreasing h implies X_i b<X_j b. If E is nonempty, the finitely many positive gaps have a positive minimum delta; divide b by delta to obtain all margins at least1. Changing the latent scale does not restrict an arbitrary response h. If E is empty, there is no constraint to satisfy. The numerical unit margin is not a measured enrichment gap or a biological effect size.

This normalization would not remain equivalent if coefficients were externally bounded or h restricted to a fixed parameterized family.

## Proof of sufficiency, including ties and unknowns

Take a feasible b and put t_i=X_i b. For each distinct latent value t, define

`y(t) = max{L_i : t_i<=t}`.

The maximum is finite. It is nondecreasing in t and at least L_j whenever t_j=t. Suppose y(t)>U_j for a row in that tied group. Then some i with t_i<=t_j has L_i>U_j. This is the forced pair (j,i), whose constraint requires t_i>=t_j+1, a contradiction. Thus y(t)<=U_j for every finite upper bound; infinite upper bounds impose no restriction. Assigning h(t)=y(t) at these finitely many scores and extending monotonically gives the required function. Piecewise-linear interpolation with constant tails can even make h continuous, though continuity is not needed.

All equal-score rows receive one common value. This covers duplicate feature rows and latent ties. Their intervals must have a common overlap, which the same forced-order test enforces. If E is empty, max_i L_i<=min_i U_i: otherwise one disjoint pair would belong to E. Then b=0 and the constant response max_i L_i work, including unknown [0,+infinity] rows. Empty input intervals must instead be rejected as invalid measurement input, not called biological null rejection.

The criterion is **not** the one for a strictly increasing h. For X=0,1,2 and intervals {0},{0},{1}, it accepts b>0 and a plateau followed by an increase. Strict increase is impossible. Conversely an arbitrary nonmonotone h would be a much broader, generally non-discriminating finite-data class. The monotonicity and shared response are substantive assumptions.

## Exact rejection certificate

Let A have rows X_j-X_i for the forced pairs. The linear alternative says A b>=1 is infeasible exactly when there is lambda>=0 with

`A^T lambda=0`, `sum(lambda)=1`.

For rational A such a certificate can be rational. Given a proposed certificate, verification is elementary and exact: every weight is nonnegative, the column sums cancel, and the weights sum to1. Multiplying the required inequalities by the weights would give 0>=1. This weighted-order cancellation rejects the whole nondecreasing codon-additive class, not merely a fitted sigmoid. Rationalizing floating weights is not itself a certificate: every identity must be rechecked, and a failed reconstruction must return unresolved. A feasible verdict likewise requires exact verification of a rational b and of the constructed response values in their intervals; optimizer status or a floating residual alone is not sufficient.

For a binary two-codon square with high00/11 and low10/01, the orders 10<00 and01<11 already cancel the contribution of the first codon: their row vectors are negatives. Weights1/2,1/2 yield an explicit contradiction. This is an elementary cancellation example and conventional reciprocal sign epistasis, not a molecular mechanism.

The test uses at most n(n-1)/2 forced pairs and standard rational linear feasibility. The intercept cancels from every difference; its presence does not add ordering freedom. The FAS factor model has8 raw additive dimensions, while these order differences have rank at most7. No special scalability claim beyond ordinary explicit rational LP is needed.

## Accessibility and interpretation

For lower/upper score fields on the *entire fixed graph*, let C_L,C_U be their endpoint widest capacities and m_L=min(L_s,L_t). Guaranteed endpoint acceptance at threshold tau is tau<=m_L. A guaranteed connection requires C_L>=tau; an optimistic disconnection requires C_U<tau. The obstruction-threshold set is exactly

`(C_U,m_L]` when C_U<m_L,

and is empty otherwise. A tied C_U=m_L gives an empty set, not a positive zero-width certificate; C_U=+infinity also gives no obstruction. An unresolved [C_L,C_U] gap is not absence or presence. Bounds here are deterministic statements about the specified envelope; they become confidence statements only if the true field has simultaneous coverage.

Empirical BR/TR envelopes are not confidence bands. Infeasibility then means that no common monotone codon-additive response fits any selection within that declared envelope. It does not establish population-level rejection, exclude technical bias, predict new genotypes or identify splicing interactions as the cause. A feasible test does not prove the true field is additive; it only fails to exclude the class. Unknown higher-order vertices can still carry outside-domain bypasses.

## Finite verification and nearest prior

`ordinal_interval_audit.py` uses exact rational Fourier-Motzkin elimination, separately enumerates all13 weak orders on three observations and checks interval compatibility directly. The two formulations share the elimination kernel, while every reported primal or dual witness is verified by direct exact algebra. Across four feature matrices and nine interval choices per row, all2916 instances agree:2646 feasible,270 infeasible, with270 algebraically verified rational duals. Cases include duplicate features, collinearity, closed ties, constant responses, finite lower/unbounded upper intervals, p=0 and empty-interval rejection. Four-point cancellation and its disappearance after interval widening are checked separately. The earlier864 finite-upper-only and v2 receipts remain historical; v3 adds an import-safe entrypoint and the explicit four-point overlap fixture. No function imports author production code. These are synthetic corroborations, not validation of a production solver or observations.

The closest original mathematical source located is Dana Scott1964, *Measurement structures and linear inequalities*, J.Math.Psychol.1(2),233-247, DOI [10.1016/0022-2496(64)90002-1](https://doi.org/10.1016/0022-2496(64)90002-1). Its publisher abstract explicitly applies finite linear-inequality solvability to measurement representations. Luce and Tukey1964, DOI [10.1016/0022-2496(64)90015-X](https://doi.org/10.1016/0022-2496(64)90015-X), is the original conjoint-measurement precedent. The full Scott source was not extractable in this bounded session: publisher403 and a reprint served as unsupported octet-stream were retained as access limitations, without claiming its exact theorem number verified. The present interval lifting and cancellation implication are supplied directly above; no novelty is inferred from access failure. Original FAS global scaling/specific-interaction work remains the relevant biological precedent, as authenticated by the separate source lane.
