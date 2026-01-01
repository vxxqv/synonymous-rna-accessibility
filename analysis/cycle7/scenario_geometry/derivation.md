# Quantifier derivation and scope

## Claim

Let G be a finite undirected graph and let f_r:V -> R_{≥0}, r=1,...,R, be a nonempty finite collection of completely specified score fields. Let P(s,t) be all simple s-to-t paths, including endpoints. Put C_r=max_P min_{v in P} f_r(v), with max over an empty path set equal to minus infinity. Put l(v)=min_r f_r(v), u(v)=max_r f_r(v), C_fix=C_l, C_ad=min_r C_r, C_max=max_r C_r and C_upper=C_u. Set E=min_r min(f_r(s),f_r(t)).

For every finite threshold tau>=0 with tau<=E, exactly one of the five categories in the protocol holds. Their exact intervals are the stated half-open threshold intervals intersected with [0,E]. The order C_fix<=C_ad<=C_max<=C_upper and the endpoint bound C_ad<=E hold. In particular, all-scenario obstruction requires C_max<tau, not C_ad<tau. This statement is proved below; the implementation is separately tested on finite synthetic ranges.

## Proof and dependencies

1. **Path criterion.** In a finite graph, connectivity in the retained vertex graph f_r>=tau is equivalent to some simple path P having min_{v in P} f_r(v)>=tau. Cycles can be removed without decreasing the minimum. Finitely many simple paths exist, so the maximum is attained when a path exists, and the capacity criterion is tau<=C_r. For an empty path set every finite threshold is disconnected. No external theorem is required.

2. **Fixed-route quantifier.** There exists one P such that, for all r and all v in P, f_r(v)>=tau if and only if there exists one P with min_{v in P} min_r f_r(v)>=tau. This is exactly tau<=C_fix. The interchange here is only between two minima, not a maximum and minimum.

3. **Adaptive-route quantifier.** For every r there exists a possibly different P_r accepted at tau if and only if tau<=C_r for every r, hence tau<=C_ad. Every fixed robust route is feasible in every scenario, so C_fix<=C_ad. This inequality can be strict. It does not assert that the route can be selected or measured without knowing which scenario applies.

4. **Some-route and no-route quantifiers.** At least one scenario connects if and only if tau<=C_max. Consequently every scenario is disconnected if and only if tau>C_max. C_ad<tau says only that at least one scenario is disconnected. It is insufficient for all-scenario obstruction.

5. **Box comparison.** Since f_r(v)<=u(v) at every vertex, bottleneck monotonicity gives C_r<=C_upper for every r. The all-scenario obstruction is therefore weaker than upper-box obstruction: C_upper<tau implies C_max<tau, but the reverse can fail. The box permits a different scenario to supply the score at each vertex, even when no joint preparation realizes that mixed field.

6. **Endpoint feasibility.** Every path contains s,t, so C_r<=min(f_r(s),f_r(t)); taking a minimum over r gives C_ad<=E. All scenarios accept both endpoints precisely when tau<=E. It is incorrect to use distinct scenario-specific endpoint thresholds while claiming a common-threshold result.

7. **Partition.** The four ordered cuts C_fix,C_ad,C_max,C_upper partition the real line into the five intervals. Their intersection with [0,E] partitions the common eligible threshold domain. Empty intervals and equal cuts cause no ambiguity. At a tie tau=C, retained vertices remain accepted, hence connectivity belongs to the interval closed on its upper boundary. If the base graph is disconnected, all capacities equal minus infinity and the upper-box obstruction category contains all of [0,E], including zero.

## Strict examples

**Adaptive routes do not imply a fixed route.** On s-a-t and s-b-t, both endpoints have score 2 in both fields. Field 1 has a=1,b=0; field 2 has a=0,b=1. Each field has capacity 1, while the vertexwise minimum has capacity 0. Thus C_fix=0<C_ad=C_max=C_upper=1 and E=2. Every threshold in (0,1] requires different routes in the two fields.

**Joint obstruction can be hidden by the box.** On the only path s-a-b-t, keep endpoints 2 and use the same two internal-score assignments. Both fields have capacity 0, while the vertexwise upper field has capacity 1. Thus C_fix=C_ad=C_max=0<C_upper=1 and E=2. All joint scenarios are disconnected at every threshold in (0,1], but the optimistic mixed box is connected.

**Mixed connectivity is not all-scenario obstruction.** On s-a-t and s-b-t with endpoints 2, field 1 has a=1,b=0 and field 2 has a=1/2,b=0. Then C_fix=C_ad=1/2 and C_max=C_upper=1. In (1/2,1], the second field is blocked and the first remains connected.

**Endpoint thresholds cannot be interchanged.** On a single edge, let both endpoints be 1/2 in field 1 and 2 in field 2. Each scenario connects at its own endpoint threshold, but a common threshold greater than 1/2 is ineligible. The common eligible domain is exactly [0,1/2].

## Unknown values

For missing entries let each unresolved score vary over [0,+infinity]. Every capacity and each of E,C_ad,C_max is monotone in every score. Thus evaluating all unresolved entries at their lower and upper endpoints gives sound metric bounds. Infinity is a symbolic upper bound, not an observed score. A certificate of one five-way category must hold for every completion. For example, adaptive-routes-only is certified by tau<=E_lower, C_fix_upper<tau and tau<=C_ad_lower. All-joint-blocked with a connected upper box is certified by tau<=E_lower, C_max_upper<tau and tau<=C_upper_lower. These are sufficient rules. Metric bounds are not assumed jointly independent, and failure to certify is unresolved, not evidence of absence.

## Non-claims and review status

This is an elementary exact accounting of existential and universal path quantifiers, not a new bottleneck algorithm or an inferred probability model. Generic minimax inequalities already explain why the two quantifiers can differ. No external prior is required to validate the internal proof, and no novelty is asserted. The distinction can be descriptive evidence only after the fixed fields are evaluated. It does not allocate a biological mechanism, identify causal epistasis, demonstrate an evolutionary route, or undo the previous global FAS exclusion or lack of an empirical-box barrier.

This derivation and its implementation are authored in one lane. The separate path-enumeration oracle is an algorithmically independent check by the same author with inherited context, not an independent external review. Root should obtain a separately authored challenge before authorizing new observed summaries.
