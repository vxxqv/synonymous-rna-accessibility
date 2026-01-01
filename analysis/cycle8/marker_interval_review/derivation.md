# Sharp marker-context range and identification bounds

## Claim card

Given a fixed finite set of M>=1 distinct, canonical marker contexts, let one unknown response y_m belong to an interval I_m=[L_m,U_m]. The response scale is declared in advance. Intervals form a Cartesian possibility set: no unrecorded cross-context constraint is assumed. Each L_m is finite or -infinity, each U_m finite or +infinity, L_m<=U_m, and each interval contains a finite real value. Infinite endpoints denote unbounded sets of finite responses, not possible infinite assay values. On a nonnegative raw-ratio scale require L_m>=0. On the log-ratio scale an unknown positive ratio is all of the real line.

Define R(y)=max_m y_m-min_m y_m. The results below describe the exact identified set conditional on the intervals, not their statistical coverage. Verdict: the closed-interval finite formulas are proved as written. Unbounded and strictly positive boundary conventions are stated separately. Proofs and fixtures here are self-reviewed pending independent challenge, not a priority claim.

## 1. The range bounds

For M=1, R is identically zero, irrespective of uncertainty in the single response. For M>=2, write

    A=max_m L_m, B=min_m U_m,
    r_min=max(0,A-B),
    r_sup=max_{m!=n}(U_m-L_n).                         (1)

There is no indeterminate infinity subtraction with the permitted endpoint conventions. In particular A cannot be +infinity and B cannot be -infinity. The upper formula uses distinct markers. For valid nonempty intervals it is nonnegative: for any two contexts, the two opposite differences sum to their nonnegative interval widths.

For bounded closed intervals the exact range identified set is [r_min,r_sup]. If at least one interval is unbounded and M>=2, r_sup=+infinity and the set is [r_min,+infinity), where every response and every realized range is finite.

### Lower bound and attainment

Every choice has max y>=A and min y<=B, so R>=A-B and R>=0. If A<=B the intervals have a common finite point, giving R=0. Select a finite point if the intersection is unbounded. If A>B, both A and B are finite. Choose y_m=max(L_m,B), interpreting max(-infinity,B)=B. It is at least L_m, at most U_m because both L_m and B are, and lies in [B,A]. A context attaining L=A gives y=A; one attaining U=B gives y=B. They must be distinct. Thus the lower bound is attained.

### Upper bound and attainment

For every choice and every distinct pair, y_m-y_n<=U_m-L_n, hence R<=r_sup. For finite r_sup, choose a maximizing distinct pair m,n and set y_m=U_m,y_n=L_n. Fill other coordinates anywhere in their intervals. This gives R>=r_sup, so equality holds. For unbounded r_sup, such a pair can be separated by an arbitrarily large finite amount. Continuity of R along line segments in the Cartesian interval domain supplies every intermediate finite range. This also proves that the bounds describe an identified interval, rather than merely two separate inequalities.

### Sharp counterexample to a naive upper bound

With two contexts I_1=[0,10] and I_2=[4,6], the exact identified range is [0,6]. The value 10 from global max U minus global min L incorrectly asks the first context to take both 0 and 10. With one context [0,10], the range is always 0, an even smaller counterexample.

### Strictly positive raw ratios and logarithms

If y>0 is required, intersect each raw interval with (0,+infinity) before use. An interval with U=0 is then empty and must fail integrity checks. For a valid finite collection the same lower bound is still attained, and (1) remains the sharp supremum. A finite upper bound need not be attained when a zero lower endpoint is excluded. For example (0,10] and (0,5] give range set [0,10), not [0,10]. If r_sup>0 is finite, it is attained precisely when some maximizing distinct endpoint pair is admissible, in particular its selected lower endpoint is positive. The r_sup=0 and M=1 cases are treated separately.

For x=log(y), use the intervals for x and compute its range anew. This is log(max y/min y), not log(max y-min y). A lower ratio bound 0 becomes an unbounded lower log bound. Numerical raw and log range thresholds are different estimands. For one wholly missing context the range remains 0; for M>=2 a missing [0,+infinity) context makes the raw and log supremum infinite. It may leave a nonzero lower bound established by the other contexts. For example [0,+infinity),[1,2],[4,5] have minimum range 2 and unbounded maximum.

## 2. A practical invariance decision

For a prospectively chosen tolerance delta>=0, define weak invariance as R<=delta. If r_sup<=delta, every completion is invariant. If r_min>delta, every completion violates invariance. Otherwise both the chosen set and its boundary conventions must be retained as unresolved; under bounded closed intervals both invariant and violating completions exist. In particular r_min=delta is not robust violation. Strict R<delta would require additional attention to whether a supremum is attained, so weak inequality is preferable for an unambiguous contract.

These are conditional sensitivity statements. Preparation envelopes do not automatically contain a population response, and empirical minima/maxima are not confidence limits. If authenticated marginal confidence sets later exist, a prespecified union bound can provide simultaneous coverage without paired covariance or independence. Neither that calibration nor a useful familywise level is supplied merely by three source-labelled rows.

## 3. What marker-only and WT controls identify

Fix one cell/amplicon and a common comparable response scale. Write

    Y_00=mu,
    Y_e0=mu+alpha,
    Y_0m=mu+beta_m,
    Y_em=mu+alpha+beta_m+eta_m.                       (2)

Here alpha is the isolated intended-edit contrast, not a name for the combined-state score. Model (2) is just a saturated contrast decomposition with baseline conventions; it assumes no absence of interaction. On a log scale alpha is an isolated log-ratio contrast. On a positive raw-ratio scale Y_e0=mu+alpha must remain positive.

Even exact WT, marker-only and combined-state observations do not identify alpha when eta_m is unrestricted:

    Y_em-Y_0m=alpha+eta_m.

For example WT=marker-only=combined=1 is compatible with edit-only response 1/2 (alpha=-1/2, eta=1/2) or 3/2 (alpha=1/2, eta=-1/2). Every observed response is positive and both directions of isolated effect remain possible. For multiple markers the simultaneous shift alpha->alpha+t, eta_m->eta_m-t preserves all those observations. On a real log scale alpha is completely unrestricted. On a positive raw scale the domain alpha>-mu remains, so it would be incorrect to claim literally every real alpha is physically admissible with fixed mu.

If C_m=[C_Lm,C_Um] bounds Y_em, D_m=[D_Lm,D_Um] bounds matched Y_0m, and a **justified absolute** interaction bound |eta_m|<=I_m is imposed, then, absent additional dependence constraints, the sharp alpha set on a real response scale is

    intersection_m [C_Lm-D_Um-I_m, C_Um-D_Lm+I_m].     (3)

Necessity follows from subtracting the controls. Sufficiency follows by selecting each difference in [C_Lm-D_Um,C_Um-D_Lm] within I_m of alpha, which is possible exactly when (3) holds; every difference in that interval is attainable by a pair of control/combined values. A shared WT interval W imposes the additional positive-raw-domain condition alpha>-W_U if its finite upper endpoint W_U is positive and attainable. On the log scale WT adds no such positivity restriction. Empty intersections mean the interval/bound assumptions are inconsistent, not a point estimate.

Formula (3) uses finite, closed, admissible intervals on the real response scale. If a raw interval's zero endpoint is excluded by strict positivity, the endpoint-admissibility conventions must also be propagated into (3); its numerical endpoints can then be infima/suprema rather than attained limits. Extended unbounded intervals are possibility sets, not finite error bars.

For known WT and only an assumed marker main-effect bound |beta_m|<=B_m, an absolute interaction bound gives the weaker necessary/sufficient per-marker alpha interval C_m-mu+[-B_m-I_m,B_m+I_m]. Both bounds are assumptions; marker-only data do not estimate I_m without edit-only information or further restrictions. Taking I_m=0 is the no-interaction model and must be labelled as such.

Marker-only controls can establish a useful relative interaction contrast without WT: for two markers a,b,

    (Y_ea-Y_0a)-(Y_eb-Y_0b)=eta_a-eta_b.              (4)

A nonzero, resolved contrast rejects equal interactions across those two marker contexts, including the zero-interaction special case. It still does not identify alpha or the absolute interaction levels. Equal adjusted contrasts can hide a shared nonzero interaction. Bounding only differences eta_a-eta_b never removes that common-shift ambiguity. Comparability of conditions/measurement is necessary, and this arithmetic is not itself a causal estimate from unpaired source rows.

## 4. Shared order across two cell contexts

For a fixed intended label and the same canonical marker states, a simple scale-robust null asks whether there are latent scores q_m and cell-specific common nondecreasing maps h_c such that y_cm=h_c(q_m). This is an ordinal marker-context statement, not an additive intended-edit model.

For exact responses in two cells, opposite strict pair orders refute this null: y_Aa<y_Ab and y_Ba>y_Bb cannot both follow one q order. Conversely, absence of such an opposite strict pair is sufficient: sort states lexicographically by (y_A,y_B), then both response sequences are nondecreasing, and define each h_c on those ordered scores. Ties are allowed and cannot be silently made strict.

For independent interval possibilities, create a directed edge a->b whenever U_ca<L_cb in either cell. There exists **some** interval completion with a common latent order if and only if this union graph is acyclic. Necessity follows because an edge forces q_a<q_b. For sufficiency take a topological order, assign strictly increasing q, and in each cell select nondecreasing responses along that order by prefix maxima of the lower bounds. A prefix lower bound exceeding the current upper bound would have forced a reverse edge, contradicting the order. For unbounded endpoints, start from a sufficiently small finite baseline rather than assigning an infinite value. This proves feasibility without fitting additive coefficients. It does not certify that the unknown true responses have that order.

Example: in cell A let the two marker responses lie in [1,2] and [3,4], and in cell B reverse those intervals. The forced two-cycle refutes common marker order on every completion. Touching intervals [1,2] and [2,3], reversed in the other cell, admit all responses equal to 2 and do not refute a nondecreasing-map null. Missing intervals can restore feasibility, so a feasible null is not evidence for invariance. Strictly increasing within-cell rescalings preserve the ordinal conclusions, while numeric range tolerances change. Arbitrary nondecreasing transformations can erase strict inequalities.

Neither reversal nor order compatibility isolates an intended edit, proves a full haplotype, identifies an RNA isoform, or resolves which bundled cell-context difference is responsible.
