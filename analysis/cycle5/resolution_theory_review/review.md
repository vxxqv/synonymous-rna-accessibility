# Deterministic resolution audit

Verdict: **proved as written**, with epsilon required nonnegative, the full Cartesian perturbation box explicitly quantified, opposite square endpoints fixed at00 and11, and weak acceptance. The statement is elementary deterministic comparison geometry. No generic theorem-priority claim is supported or needed.

## Claim card

The graph is the undirected four-cycle with vertices00,01,10,11 and single-coordinate edges. Its endpoints are s=00 and t=11. The two trait vectors are f(s)=f(t)=(0,0), f(01)=(m,-m) and f(10)=(-m,m), where m>0. For each vertex v and trait k, choose any real error e(v,k) in the closed interval[-epsilon,epsilon], with epsilon>=0. All eight errors range jointly over the full Cartesian box; this is a deterministic set, not a probability or independence assumption. Put f'=f+e and recompute tau'_k=min(f'_k(s),f'_k(t)). Acceptance is f'_k(v)>=tau'_k, including equality.

The claim is that every such perturbation is joint-only obstructed, meaning each individual trait connects s to t but the intersection does not, if and only if m>2epsilon. Endpoints are automatically jointly acceptable by this definition.

## Proof and sharpness

For m>2epsilon, each threshold lies in[-epsilon,epsilon]. The originally positive intermediate coordinate is at least m-epsilon>epsilon>=tau'_k, so each trait retains its two-edge route. The originally negative coordinate is at most -m+epsilon<-epsilon<=tau'_k, so each intermediate fails one of the two traits. Every s-to-t path on this square must contain at least one of its two intermediate vertices. Neither intermediate is jointly acceptable, and hence no joint path exists. Repeated vertices cannot bypass these two vertices. This argument covers arbitrary real errors in the box, including all interiors; no corner enumeration is needed for sufficiency.

For 0<m<=2epsilon, set both endpoints to(-epsilon,-epsilon). Perturb01 to(m-epsilon,-m+epsilon), using allowed errors(-epsilon,+epsilon). At the recomputed threshold(-epsilon,-epsilon), its first coordinate passes because m>=0, and its second passes exactly when m<=2epsilon. Thus s,01,t is a joint route. The other intermediate can be left unchanged. This is an admissible box perturbation, so the universal assertion fails. At equality its second coordinate equals the threshold; the weak acceptance rule is decisive. At epsilon=0, all m>0 satisfy the assertion. The excluded m=0 field has an unperturbed joint route.

Swapping s and t or relabeling the whole square together with its field preserves the result. Choosing the other diagonal in the same fixed field does not: its endpoint minimum is(-m,-m), so both(0,0) intermediates pass and the unperturbed diagonal is jointly connected. Adjacent jointly acceptable endpoints are connected by their edge. The opposite-endpoint hypothesis cannot be omitted.

## Comparison perturbation bound

For any non-endpoint v and trait k, write g=f_k(v)-min(f_k(s),f_k(t)). If all three errors have absolute value at most epsilon_k, then

```text
|min(f'_k(s),f'_k(t))-min(f_k(s),f_k(t))| <= epsilon_k,
|g'-g| <= 2epsilon_k.
```

The first inequality follows by bounding each perturbed endpoint above and below by its original value plus or minus epsilon_k and then taking minima. The second is the triangle inequality. The constant2 is attained by moving v upward and both endpoints downward, or vice versa, so it cannot be reduced for arbitrary Cartesian errors. Endpoints themselves are always acceptable; treating their acceptance as an independent sign comparison would be unnecessary and could be misleading.

A sufficient way to preserve a given finite acceptability pattern is to keep every rejected comparison strictly below-2epsilon_k and every accepted comparison at least2epsilon_k above zero. The more conservative formulation that all signed margins are strictly greater than2epsilon_k is also sufficient. It is not necessary. For example, with epsilon=1, endpoints(0,0), intermediates(2,-3) and(-3,2), every perturbation still has the intended positive coordinates weakly acceptable and the negative coordinates strictly rejected. The positive margin equals2epsilon rather than exceeding it. Comparisons at irrelevant vertices can also have small margins without controlling a particular connectivity conclusion. These are sign-pattern sufficient conditions, not a general necessary-and-sufficient graph characterization.

## Fixed-low-end interval certificate is a different question

For the symmetric field, give every coordinate its interval[f-epsilon,f+epsilon]. The protocol's separate fixed threshold is tau_k=min(L_k(s),L_k(t))=-epsilon. Each individual pessimistic graph has its intended positive intermediate, since m-epsilon>=-epsilon. Its negative coordinate's upper bound passes the optimistic joint test exactly when -m+epsilon>=-epsilon, or m<=2epsilon. Consequently **this symmetric example's fixed-low-end range certificate also passes exactly when m>2epsilon**. It would be incorrect to claim a different threshold for this example.

The apparent equality does not identify the two definitions in general. A decisive counterexample has epsilon=1, endpoints(0,0), and intermediates(1,-3),(-3,1). At fixed tau=(-1,-1), the individual lower-bound graphs each have a route, while the negative-coordinate upper bounds are-2<-1, so the optimistic joint graph is disconnected. The range certificate passes. But the admissible perturbation

```text
endpoints: (1,1), (1,1)
intermediates: (0,-4), (-4,0)
```

has endpoint-relative threshold(1,1), and both individual trait graphs are disconnected. Therefore the fixed-low-end certificate does not certify universal joint-only obstruction when endpoint-relative thresholds are recomputed for every realization. Its correct scope is robustness at its separately fixed threshold.

For independent coordinate intervals on this square, a universal endpoint-relative joint-only obstruction does imply the fixed-low-end certificate: set endpoints to their upper bounds and both intermediates to their lower bounds to obtain the individual pessimistic routes, which also pass at the lower fixed threshold. Set endpoints to their lower bounds and both intermediates to their upper bounds to see optimistic joint disconnection. These are two distinct admissible realizations. The reverse implication fails by the example above. This observation is confined to the declared coordinate-interval interpretation, and is not a new claim about measured uncertainty coverage.

## Common offsets, batches and interpretation

Adding the same c_k to trait k at every vertex adds c_k to its endpoint minimum. Every comparison is unchanged exactly. More generally e(v,k)=c_k+d(v,k) removes c_k algebraically; a deterministic residual bound |d(v,k)|<=eta_k yields the comparison bound2eta_k. A batch or family error shared by only some vertices need not cancel. Endpoints in one batch and intermediates in another can move their relative comparisons. A covariance description alone neither supplies a deterministic error box nor establishes that the relevant errors are a common offset. A restricted correlated uncertainty set can omit the adversarial box witnesses, so the full-box necessity direction must not be transferred to that smaller set. For pure shared offsets, every m>0 preserves the symmetric obstruction, however large the offset bound.

No measured calibration of epsilon is assumed here. Empirical preparation ranges do not automatically contain true latent scores, and their twelve dependent combinations do not create a confidence level. The numerical tie diagnostic is not a noise bound. This deterministic result therefore supplies neither a confidence interval nor empirical power, a noise ceiling, or exclusion of small true effects when observed obstructions are absent. Its useful implication is conditional and explicit: a symmetric conflict of signed size m is immune to a declared independent bounded-error set exactly below the radius m/2, with strict radius because weak boundary ties can create a route.

## Dependencies and audit obligations

The dependency graph consists only of the square adjacency convention, weak acceptance, a minimum's elementary uniform bound, the full Cartesian error set, and the explicit witnesses above. There are no unchecked external theorem leaves. Definitions, arbitrary-interior sufficiency, sharpness, equality, epsilon=0, endpoint orientation, interval distinction and shared-offset cancellation pass. Biological calibration and general full-graph characterization are out of scope.

`verify.py` independently implements graph traversal using exact rational arithmetic. It checks a declared ternary error grid and explicit witnesses. The grid is finite evidence only: the continuum conclusion rests on the proof above, not on grid coverage. The script refuses optimization and overwrite. The command and before/after source hashes are retained in `run/receipt.json`. No biological files, graph implementation from the pilot, measured outcomes, manuscript or repository files are read. This is an independent challenge of the parent's proposed claim, not a review of an empirical finding.
