# Finite sensitivity and missing-label inference

This is an exploratory extension of the exposed B1 audit. The proofs below use standard finite order, interval and enumeration arguments. The application-specific results require separate numerical review before manuscript admission.

## Exact finite-marker transition

For M>=2 strictly positive finite closed intervals I_m=[L_m,U_m], let A=max_m L_m/min_m U_m and B=max_{m!=n}U_m/L_n. Under the declared common widening I_m(gamma)=[L_m/gamma,gamma U_m], t=gamma squared>=1, the established sharp dispersion formulas give

F_min(t)=max(1,A/t), F_max(t)=Bt.

Equality is compatible iff t>=max(1,A), including the boundary. Minimum dispersion at least2 holds iff1<=t<=A/2, provided A>=2. Maximum dispersion at most2 holds iff1<=t<=2/B, provided B<=2. The strict/equality distinction is essential. These statements follow by substitution, not by a new generic optimization theorem. Unknown intervals require the original extended endpoint conventions; the current complete source has none, so this follow-up does not assert finite formulas for unknown endpoints.

## Exact unrestricted common-order transition

Let V be the finite set of matched local states in an exon and c its contexts. In scaled intervals J_cv(t)=[L_cv,t U_cv], a common total ordering is feasible iff the graph containing i->j whenever some context has t U_ci<L_cj is acyclic. Necessity follows because separate nondecreasing responses of a shared scalar cannot reverse a strict forced order. For sufficiency take any topological ordering v_1,...,v_n and, for each context, assign y_c(v_j)=max_{i<=j}L_c(v_i). If this exceeded t U_c(v_j), one earlier i would force v_j->v_i, contradicting the topological order. Thus all responses lie in their scaled intervals and are nondecreasing. Giving the ordered states distinct integer latent scores and interpolating by steps constructs the required response functions. Division by gamma converts to the original widened intervals. This construction permits ties in responses and uses strictly positive lower endpoints.

Write w_ij=max_c L_cj/U_ci. A directed cycle remains present precisely when t<min_{edges e in cycle}w_e. Consequently

t_star=max(1,max_directed_cycles min_{edges e in cycle}w_e),

where the empty cycle set contributes no value above1. This finite maximum is attained. The graph is acyclic at t_star because edges with w_e=t_star are absent at equality. If t_star>1, a cycle whose minimum edge weight is t_star certifies infeasibility for every1<=t<t_star. A topological order and rational scaled responses certify feasibility at t_star. Searching the finitely many rational weights is exact and terminates. No floating-point root or tolerance decides feasibility.

This concerns an arbitrary shared scalar assigned to each state. It neither proves nor fits an additive sequence score. Acyclicity cannot be used to claim additive representability. The maximum-cycle-bottleneck argument is elementary graph reasoning; the distinct contribution is the exact boundary for these specified measurements and this specified uncertainty family. It is not an empirical estimate of measurement error.

## Anonymous preparation assignments

Let each of M measured marker states supply three paired RNA/DNA ratios. The values within each original row are inseparable. Unknown cross-marker labels permit every bijection of the three rows of each marker to three common labels. Fixing the first marker's labeling removes only a simultaneous relabeling symmetry: every full assignment is equivalent to exactly one represented assignment for distinct row identities. There are(3!)^(M-1) representatives, even when some numerical values coincide, because row identities remain distinct. Summaries invariant to simultaneous relabeling therefore have exactly the same possible range under this enumeration as under all(3!)^M labelings.

For each representative compute three within-label dispersions, count those at least2, and take the minimum and maximum of this count over representatives. The lower count holds for every admissible local assignment. The upper count is attainable for at least one. They are finite recorded-data facts under the weak assignment set, not confidence bounds on true preparation effects. They do not authenticate the labels or independence of cultures. Cross-group source constraints could narrow the assignment family if recovered.

For two markers, a separate combinatorial characterization is available. Make a bipartite graph of the three rows on each side, marking an edge good when its ratio dispersion is less than2. The minimum number of at-least-twofold matched pairs is3 minus the maximum matching size in the good-edge graph. A maximum good matching can always be completed to a full bijection in the complete bipartite graph. Conversely a bijection cannot contain more good edges than a maximum matching. The maximum at-least-twofold count is the maximum matching size in the complementary bad-edge graph. This provides an independent checker of all two-marker cases. It is a standard matching consequence.

## Important non-implication

Point-set assignment robustness is different from continuous-interval uncertainty robustness. For invented marker triples(1,4,16) and(2,8,32), every recorded cross-marker pair differs by at least2, so every assignment yields three at-least-twofold comparisons. Their continuous intervals[1,16] and[2,32] nevertheless overlap, so finite-envelope F_min=1. The extra recorded-point statement cannot be advertised as stronger biological confidence or as a correction of the primary envelope result. Opposite directions across preparation labels may also contribute to dispersion and do not identify a consistent signed marker effect.

## Closest-prior boundary

Robust ordinal regression already distinguishes necessary from possible statements over compatible models and uncertain evaluations. Interval order extension, cycle detection and matching are established tools. The present derivations adapt those tools to exact deposited ratios and explicit lost-label ambiguity. No claim is made that a source's marker non-neutrality, a generic rank method, or these elementary graph identities are newly discovered. The biological evidence remains a selected local-motif RNA/DNA signal from one study.
