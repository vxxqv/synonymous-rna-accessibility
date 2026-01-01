# Exact full-system sensitivity thresholds for an ordinal additive representation

Primary verdict: **correct only after stated restrictions and certificate qualifications**. The proposed finite-breakpoint characterization is correct for the rational finite-domain common-nondecreasing-response model. Reference anchoring changes the breakpoints. Zero upper bounds or incompatible fixed constraints can give infinite thresholds. An arbitrary infeasible certificate somewhere below a proposed boundary does not by itself prove that boundary is sharp.

This is a direct sensitivity corollary of the Cycle6 interval/order characterization, not a foundational ordinal-representation or robustness novelty claim. It does not calculate any observed successor threshold.

## Claim card

Fix n>=1 rational feature rows X_v in Q^p, with p>=0. The feature design is fixed independently of the sensitivity parameter. Let 0<=L_v<=U_v with finite rational L_v; U_v may be rational or +infinity. Assign a fixed binary expansion indicator a_v in {0,1}. A reference anchored at its normalized value1 has L_w=U_w=1 and a_w=0. Other vertices have a_v=1. Multiple fixed vertices are allowed in the general statement but can cause permanent infeasibility.

At each parameter, ask whether there exist beta in R^p and one common nondecreasing real function h such that h(X_v beta) belongs to every expanded interval. The function is not required to be strictly increasing, sigmoid, differentiable, or bounded in slope. All intervals constrain nonnegative outputs. The beta and h may change with the parameter. The question is existential representability of the full declared domain, not preservation of one fitted model.

Additive expansion at epsilon>=0:

`I_v(epsilon)=[max(0,L_v-a_v*epsilon), U_v+a_v*epsilon]`.

Symmetric multiplicative expansion at gamma>=1:

`J_v(gamma)=[L_v/gamma^a_v, U_v*gamma^a_v]`.

Infinite upper bounds remain infinite. The threshold is the smallest feasible parameter if one exists; otherwise it is +infinity. The main scientific metric fixes the normalized reference. The all-a=1 expansion is a different, unanchored secondary metric.

## Dependency graph and evidence status

1. Definitions and nonnegative, closed intervals: conventions and explicit assumptions, passed.
2. Finite interval/order equivalence, proved in the reviewed Cycle6 characterization: internal prior dependency, reconstructed below, passed.
3. Exact order-expiry formulas and nested active sets: direct algebra, passed.
4. First feasible breakpoint and attainment: finite-set argument, passed.
5. Exact feasible and dual verification: direct algebra plus the classical rational linear-inequality alternative, passed as a finite certificate characterization. Floating solver success alone is not evidence.
6. Synthetic checks:3136 exact comparisons with a separate rank-one weak-order oracle;866 squared-interval equivalence checks. These are limited implementation checks, not the proof.
7. Biological calibration, sampling coverage and absence of persistent recovery bias: not established and out of scope.

## Interval/order equivalence

For any valid intervals with finite lower bounds, representation is equivalent to the rational linear system

`(X_j-X_i) beta >=1` for every strict forced order `U_i<L_j`.

Necessity: output separation requires strictly ordered latent scores; finitely many positive latent gaps can be rescaled to have minimum1. No forced orders require no rescaling. Sufficiency: for latent scores z_v from a feasible beta, assign `y_v=max{L_j:z_j<=z_v}`. This is nondecreasing and constant at ties. If y_i>U_i, some preceding lower bound L_j>U_i would force z_i<z_j, a contradiction. Infinite upper bounds cannot violate the upper check. Extend this finite assignment by monotone interpolation and constant tails. With no forced orders, all intervals share max L_v and a constant response works. The argument also covers p=0 and duplicate feature rows.

## Additive order expiries and sharp threshold

Put k_ij=a_i+a_j. Since U_i>=0,

`max(0,L_j-a_j*epsilon)>U_i+a_i*epsilon`

holds exactly when `L_j-U_i>k_ij*epsilon`. If the clipped lower bound is0 it cannot exceed the nonnegative upper bound; otherwise clipping does not change the comparison. No outgoing order is contributed by U_i=+infinity.

For k_ij>0, an initially forced order disappears precisely at

`d_ij=(L_j-U_i)/k_ij >0`.

Equality is not a strict order. For k_ij=0, an initially forced order persists forever. Thus active systems are nested and constant between the finite sorted candidates `{0} union {d_ij>0}`. The smallest feasible epsilon, if finite, is one of these rational candidates, and it is attained. Indeed, the system at an event already omits all orders expiring there, and stays unchanged until the next event. All equal-time removals must occur together. If the final permanent-order system is infeasible, no finite expansion works; otherwise the final candidate works.

For the all-a=1 metric, d_ij=(L_j-U_i)/2. With a single anchor[1,1] and every other vertex expanded, a sufficiently large finite epsilon makes every interval contain1, so a constant response proves finite existence. Without an anchor the same is true for sufficiently large epsilon and response0. No assumption about graph connectivity or epistasis enters this algebra.

## Multiplicative order expiries and exact construction

Write t=gamma^2>=1, not gamma. A forced order is

`L_j > U_i*gamma^(a_i+a_j)`.

For U_i>0 and L_j>U_i, its finite removal candidate is:

- k_ij=2: `t_ij=L_j/U_i`;
- k_ij=1: `t_ij=(L_j/U_i)^2`;
- k_ij=0: no removal, hence permanent.

If U_i=0 and L_j>0, the order is also permanent. If L_j=0 it cannot be an outgoing forced comparison from a nonnegative U_i. Unknown infinite upper bounds give no outgoing comparison. Every finite t candidate is rational. The same nested finite-set argument proves that the first feasible t is an attained candidate in `{1} union {t_ij>1}`, unless the permanent system is infeasible. Then the threshold is infinite. The smallest gamma is sqrt(t_star), which need not be rational.

An exact rational boundary construction is available even for non-square t_star. Square the physical expanded intervals:

`K_v(t)=[L_v^2/t^a_v, U_v^2*t^a_v]`.

All these endpoints are rational when t is rational, and the reference remains[1,1]. Squaring preserves nonnegative order, so the original representation exists iff one common nondecreasing function g maps the latent scores into K_v(t). For one direction, replace h by max(h,0) if necessary and put g=max(h,0)^2; the observed nonnegative values are unchanged. For the reverse direction put h=sqrt(g), using a nonnegative monotone extension of the finite g values. A rational beta and rational prefix-max responses g_v are therefore an exact feasible certificate. Physical responses are the explicitly specified roots sqrt(g_v); no irrational optimizer or rounded gamma is required.

Squaring does not change feature-difference rows, forced-order identities or dual weights. In particular it does not change the scientific metric from gamma to a different expansion. It is only an algebraic representation of the same interval membership. The nonnegative assumption is essential: squaring arbitrary signed outputs would not preserve order.

For positive finite U_v, the unanchored threshold is finite. With one reference1 fixed and all other vertices expanded it is also finite: sufficiently large gamma makes every nonreference interval contain1. Multiple fixed vertices or zero upper bounds require the permanent-system check and can produce infinity.

## What proves a proposed threshold sharp

Upper certificate: at candidate epsilon_star or t_star, verify a rational beta against every active order and construct the monotone interval responses. For the multiplicative boundary use the squared rational intervals and record response_squared, not rounded physical outputs. Include every vertex and all forced pairs.

Lower certificate: supply nonnegative rational weights lambda with sum1 and exact feature cancellation `sum lambda_ij*(X_j-X_i)=0`. Each positively weighted comparison must remain strict for **every smaller admissible parameter**. Equivalently, its expiry must be at least the proposed boundary, or it must be permanent. Then their weighted latent margins would sum to a positive number at any smaller parameter although the feature sum is0. This excludes every smaller parameter without rerunning all intermediate LPs.

For an additive expiring order the lifetime condition is `L_j-U_i >= k_ij*epsilon_star`. For multiplicative k=1 it is `L_j^2>=U_i^2*t_star`; for k=2 it is `L_j>=U_i*t_star`. Permanent orders must still be strictly forced. Check the comparison's type and positivity, rather than accepting a non-strict inequality for a pair that was never forced.

A dual obtained in the immediately preceding breakpoint cell automatically has these lifetimes. An arbitrary dual earlier in the search need not. If the initial parameter0 or1 is already feasible, its threshold is the initial parameter and no lower-domain certificate is needed. If no finite parameter is feasible, a dual of the final permanent system certifies infinity.

Counterexample to an insufficient sharpness claim: two identical feature rows with point responses1 and3 need epsilon=1. A dual at epsilon0 and a feasible construction at epsilon2 are both valid, but do not make2 sharp. Also, loss of the first selected witness need not restore the whole system: duplicate pairs with point responses(1,2) and(1,4), at feature rows(0,0,1,1), have first witness expiry epsilon1/2 but full-system threshold3/2; corresponding unanchored multiplicative t expiries are2 and4. These are exact synthetic values, not FAS measurements.

## Boundary controls and known answers

The file synthetic_receipt_v2.json records these invented controls. X=(0,1,2) is a one-dimensional latent design:

- Point outputs(1,4,2), all expanded: epsilon_star=1 and t_star=2. At the boundary the latter two outputs can tie. gamma_star=sqrt(2), so an all-rational-gamma requirement would wrongly omit the exact optimum.
- Same points with the middle vertex fixed: epsilon_star=2 and t_star=4. Anchoring changes the answer.
- All outputs1: epsilon_star=0, t_star=1, constant response.
- Point outputs(0,1,0), all expanded: epsilon_star=1/2 but no finite multiplicative factor. The two zero-output endpoints stay zero, so a monotone response on three ordered latent positions cannot have a positive middle peak for either latent direction.
- All intervals[0,infinity]: initial constant response is feasible.
- Two identical feature rows fixed at different values: neither additive nor multiplicative widening of other vertices repairs their permanent contradiction.

Changing nondecreasing to strictly increasing would invalidate plateau constructions. Negative or empty intervals are rejected as model-integrity failures, not silently clipped into a different initial problem. If strict error budgets epsilon<e rather than <=e are intended, boundary attainment is lost and the decision must be stated separately.

## Complexity and interpretation

There are at most n(n-1) directed candidate orders, O(n^2) event values, and a rational LP in p variables per distinct active set. Exact sorting and binary search can reduce decision calls to O(log n), plus certificate extraction and verification. This is a polynomial-bit rational-LP statement with the input bit length included, not a strongly polynomial guarantee or a runtime bound for the toy Fourier-elimination code. That code is deliberately restricted to small synthetic cases. Irrational input intervals or features lose the stated rational-certificate guarantee; general nonbinary expansion powers can also make t breakpoints irrational.

Below the sharp threshold, every field in the stated expanded box remains incompatible with the common monotone additive model. At and above it, at least one compatible field exists. This is not the claim that every perturbation repairs the model, that the actual unknown error has the required size, or that a fitted h identifies any causal process. A given tuple of errors may be constrained by source dependence; the Cartesian relaxation can admit combinations not observed jointly. Shared global scale changes and genotype-specific recovery biases are different, and an anchored metric must not be replaced by an unanchored one without changing its label.

A positive radius measures distance from the declared empirical interval set to the finite-domain model class under the chosen metric. It is not a confidence level, measurement calibration, interaction strength, or demonstrated robust biological mechanism. Persistent genotype-dependent recovery can move the true endpoint outside every proposed empirical envelope; the radius alone cannot distinguish such bias from nonadditivity in the measured biological process. No graph obstruction or function improvement follows from this calculation.

## Review status and remaining step

This derivation is a new internal proof supplied in response to the parent's requested audit, with synthetic self-checks and no observed successor calculation. It needs the separate fresh adversarial review requested by the parent before manuscript use. The cheapest next check is to freeze the anchored metric, all event/tie rules and squared-certificate schema, then have an independent verifier inspect both sharpness certificates; no new observed subset or threshold selection is justified by this proof.
