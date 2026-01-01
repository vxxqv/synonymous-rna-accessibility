# What two decay measurements identify

**Primary verdict: point identification of a pathway contribution from two half-lives alone is refuted.** The strongest directly identified quantity, conditional on valid comparable single-exponential half-lives, is the net change in total decay rate. Residual activity and background shifts require explicit restrictions before that change can be assigned to a pathway. These are elementary competing-rate consequences, not a novel kinetic theorem or a biological result.

## Contract and measured quantities

Assume `N(t)=N(0) exp(-kt)` with a positive, time-invariant rate, and half-lives `t0,t1` on the same identified time scale. Then

`k0=ln(2)/t0`, `k1=ln(2)/t1`, `D=k0-k1`.

This transformation is conditional on the endpoint being a half-life. A steady-state abundance, unspecified decay coefficient, mixed-exponential fit or differently normalized number cannot simply be substituted. No measured calibration, half-life fit or CNOT3 outcome is assumed or examined here.

Consider

`k0=b+p`, `k1=b+delta+r*p`, `b,p>=0`, `0<=r<=1`.

Here `p` is the baseline pathway-associated rate, `r*p` its residual contribution, and `delta` the signed change in the aggregate other/background rate. Calling r a residual fraction is an assumption; r>1 would describe activation and falls outside the bounds below. Positive delta is compensatory acceleration of the other channel; a negative delta is a background slowdown. Define the suppressed pathway rate `S=(1-r)p`. The observable identity is

`D=S-delta=(1-r)p-delta`.                                    (1)

The model as initially written requires only b,p>=0. If b+delta is itself interpreted as a physical background rate, add **b+delta>=0**. That extra condition is stated separately below and in code. Every nonidentifiability counterexample here satisfies the stronger physical condition too.

## Exact rates: identification, boundaries and contradictions

If r<1 and delta are known, then

`p=(D+delta)/(1-r)`, `b=k0-p`.                                (2)

This is a feasible point identification only if `0<=p<=k0`. The physical version additionally requires `r*p<=k1`. A negative inferred p or an upper-bound violation rejects the joint restrictions at those exact rates; it must not be clamped into a different model.

If r=1, equation (1) requires `delta=-D` and contains no p. When that equality holds, the bare model permits every p in [0,k0]; the physical version permits [0,min(k0,k1)]. Boundary intervals can collapse to a point, but r=1 does not generally identify a pathway rate.

With no background shift, delta=0, and known r<1:

- D>0 implies p>0, but p=D only for complete suppression r=0.
- D=0 implies p=0. If r=1 is also admissible, D=0 can instead hide any contribution compatible with the baseline.
- D<0 is incompatible with this no-shift inhibitory model. It is not evidence for a negative pathway rate. Signed background shifts, residual-activity assumptions, endpoint comparability or measurement error remain competing explanations.

With allowed compensation, zero response can hide a positive contribution: b=1,p=2,r=1/2,delta=1 gives k0=k1=3. A rate increase can coexist with that same p: delta=2 instead gives k1=4. Conversely b=2,p=0,r=0,delta=-1 gives k0=2,k1=1: the positive net reduction comes entirely from the background shift.

Two decisive point-identification failures are:

| Shared exact observables | First (b,p,r,delta) | Second (b,p,r,delta) |
|---|---|---|
| k0=2, k1=1; delta=0 | (1,1,0,0) | (0,2,1/2,0) |
| k0=3, k1=2; r=1/2 | (1,2,1/2,0) | (0,3,1/2,1/2) |

The first varies residual activity; the second holds it fixed and varies only nonnegative compensation. Both reproduce the same measured rates with different p.

## Sharp exact-rate bounds with externally stipulated uncertainty

Let r lie in [rL,rU] subset [0,1] and delta in [dL,dU], with every combination allowed. Write `qL=1-rU`, `qU=1-rL`. The complete bare-model identification set consists exactly of p satisfying

`0<=p<=k0`, `qU*p>=D+dL`, `qL*p<=D+dU`.                       (3)

The physical version adds

`rL*p<=k1`, `p<=k0+dU`.                                      (4)

These are scalar linear inequalities, including r=1 cases. A zero coefficient is a feasibility condition, not an invitation to divide by zero. For example, `0*p<=D+dU` requires D+dU>=0. Empty sets are retained as model incompatibility.

To prove sharpness, put w=(1-r)p=D+delta. For a fixed candidate p, w must belong to both [qL*p,qU*p] and [D+dL,D+dU]. Their intersection is nonempty exactly under (3). The physical constraint adds w>=p-k1; intersection with that half-line gives precisely (4). For p>0 choose any allowed w, set r=1-w/p, delta=w-D and b=k0-p. For p=0 choose delta=-D and any allowed r. This constructs a witness for every retained p, not just the extrema.

For fixed r<1 and bounded delta, the bare interval simplifies to

`[max(0,(D+dL)/(1-r)), min(k0,(D+dU)/(1-r))]`,

provided the lower endpoint does not exceed the upper. The physical version also caps p at k1/r when r>0. For delta=0 and positive D with completely unknown r in [0,1], the sharp interval is [D,k0], and the observations require `r<=k1/k0`. Thus incomplete inhibition makes D a lower bound, not generally the baseline pathway rate. Nonnegative compensation preserves the inequality D<=S<=p; an unconstrained negative delta removes it.

## Deterministic rate intervals

Suppose calibrated rates belong to finite intervals [L0,U0] and [L1,U1], with nonnegative endpoints, and r/delta have the independent bounds above. The sharp p set for this declared rectangular admissible set is exactly

`0<=p<=U0`,

`(1-rL)*p >= L0+dL-U1`,

`(1-rU)*p <= U0+dU-L1`,

`rL*p <= U1-dL`.                                             (5)

The physical version adds `p<=U0+dU` and `rL*p<=U1`.

For proof, fix p. Its possible b values are [max(0,L0-p),U0-p]. Let n=b+delta. In the bare model n fills the interval [max(0,L0-p)+dL,U0-p+dU]; in the physical model intersect this interval with [0,infinity). Independently, c=r*p fills [rL*p,rU*p]. Feasibility is exactly intersection of the resulting n+c interval with [L1,U1]. Expanding its endpoint inequalities gives (5) and the two physical additions. Each interval-sum step is attainable. The reusable code explicitly constructs b,delta,r and both rates for any retained p, supplying a sharpness witness.

For positive half-life intervals `t_j in [a_j,z_j]`, rate intervals reverse their endpoints: `[ln(2)/z_j,ln(2)/a_j]`. Infinite upper half-lives have limiting rate lower bound zero; if only finite lifetimes are admitted, that endpoint is an infimum rather than an attained rate. A zero or negative half-life is not an admissible finite input. A band crossing a contradictory sign may restore feasibility; an exact contradictory estimate is never silently replaced by zero.

These are deterministic identified sets. They are not confidence intervals unless an independently justified simultaneous measurement statement supplies that interpretation. Correlation among rate, residual-activity or compensation bounds can make the admissible set smaller; the rectangular formulas then remain outer bounds but need not be sharp. A common unknown positive time-scale multiplier preserves fold changes and signs/rankings of absolute rate differences but leaves their absolute units unidentified; condition-specific unknown scale factors can invalidate even cross-condition comparisons.

## Fold changes and baseline coupling

For finite positive rates, the half-life fold change is

`H=t1/t0=k0/k1`, `D=k0*(1-1/H)`.                             (6)

H>1 is equivalent to D>0, but H alone supplies no absolute rate reduction. The rate pairs (1,1/2) and (10,5) both have H=2 and reductions 1/2 and 5. Rankings can reverse: (4,2) has H=2,D=2, while (1,1/4) has H=4,D=3/4.

Under delta=0 and known r<1, the baseline pathway fraction is `p/k0=(1-1/H)/(1-r)`. In particular H<=1/r when r>0; a larger fold is incompatible with that residual fraction and the no-shift model. With unknown r, a fold change is not a point estimate of that fraction.

The baseline appears algebraically in the response: with constant pathway fraction s=p/k0 and common r,delta=0, `D=(1-r)*s*k0` even though fractional dependence is identical. With constant p and r but varying b, D is constant while `H=k0/(k0-D)` grows as the baseline rate decreases when D>0; at D=0 the fold is constantly one. Baseline-associated features can therefore correlate with a response through its definition or varying background. A shared error in the control rate enters both the reported baseline and D with the same sign. No error distribution or empirical calibration is inferred here.

**Useful CNOT3 interpretation consequence:** stable baseline RNA and a positive pathway-associated decay contribution are not contradictory. For example b=1/100,p=9/100 gives total rate 1/10, lower than a different RNA with b=1,p=0. Cross-transcript background variation can reverse baseline stability even when the pathway increases decay holding b fixed. Therefore a stable-transcript enrichment does not by itself refute pathway-associated decay, and a large half-life fold response does not by itself rank the absolute pathway contribution. A defensible analysis must distinguish baseline lifetime, net rate reduction, suppressed contribution and baseline pathway contribution, while keeping r/delta assumptions visible. This is an interpretive consequence to test on eligible measurements, not a result about the CNOT3 data.

## Bounded verification, dependency limits and stopping point

`python -B verify.py --out NEW_DIRECTORY` runs 24 hand-declared rational fixtures under both model readings: 48 cases, 12,960 exact four-constraint systems and 94 constructive witnesses. Every hand-derived interval agreed with an independently implemented four-variable polytope-vertex oracle based directly on k0=b+p and k1=b+delta+c. Four observationally identical/different-p counterexamples, same-fold/different-drop cases and a fold/drop rank reversal were checked exactly. No random, floating or biological input was used. The passed receipt is `synthetic_v1/results.json`. These finite checks corroborate the algebra in this stated range; they are not exhaustive parameter coverage or the proof of sharpness.

Dependency chain: declared single-exponential/comparability assumptions -> reciprocal transformation -> identity (1) -> explicit parameter restrictions -> identification/bounds -> endpoint interpretation. Algebra, boundary handling, witness construction and named synthetic checks pass. Biological exponential adequacy, rate calibration, CNOT3-specific r/delta restrictions, measurement uncertainty and causal decomposition remain conditional/unexamined. Prior attribution is conservative: competing rates, reciprocal half-life conversion and partial identification by elementary elimination are established tools; no new generic theorem or empirical effect is claimed.

All work is confined to this model directory. No biological outcome, shared manuscript or repository was read for numerical values or modified; no commit, publication, experimental design or engineering action occurred. This derivation and the bounded counterexamples complete the assigned static inference task.
