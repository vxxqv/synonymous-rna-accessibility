# Positive steady-state identification fiber

## Claim card and verdict

**Verdict: proved as written for the mathematical steady-state observation model below.** This is identification from one exact equilibrium, not a claim about identifiability from complete time trajectories. The derivation is newly authored here and self-checked; it has not yet received the fresh independent review required before manuscript use.

Let s,e,d_N,d_C,k,d_P be strictly positive, finite, constant parameters in

    N' = s - (e+d_N)N,
    C' = eN - d_C C,
    P' = kC - d_P P.

N and C are calibrated, compatible amounts of the same nuclear/cytosolic RNA population. P is calibrated amount of the protein product. Parameters are effective rates within this model; k is completed protein production per cytosolic RNA per time, not automatically an initiation rate. No reverse import, feedback, maturation delay, cell-state mixture or separate translation-inactive RNA pool is represented. Dilution, if relevant, must be consistently included in the effective losses.

For every observed triple (n,c,p) in the strictly positive orthant, the **entire** parameter fiber is

    e > 0, d_N > 0, d_P > 0 arbitrary,
    s = n(e+d_N),
    d_C = en/c,
    k = d_P p/c.                                      (1)

### Necessity, sufficiency and stability

At equilibrium each derivative is zero, so s=(e+d_N)n, en=d_C c and kc=d_P p. Division is legitimate because n,c,p are positive. These equations force (1) and do not further restrict its three free parameters. Conversely every positive choice in (1) sets all derivatives to zero. The affine system's linear part is lower triangular with eigenvalues -(e+d_N), -d_C and -d_P. All have negative real part, including repeated-eigenvalue cases, so this equilibrium is unique and globally attracting for this model. This dynamical fact does not establish that any particular assay had reached equilibrium.

The map from parameters to (n,c,p,e,d_N,d_P) and its inverse (1) are smooth on their positive domains. Consequently the fiber has exactly three free continuous dimensions. This is stronger and more precise than merely exhibiting one confounded pair.

## What is, and is not, identified

An equivalent parameterization uses a>0, f in (0,1), and b>0:

    s=na, e=fa, d_N=(1-f)a, d_C=fan/c,
    d_P=b, k=bp/c.                                    (2)

Here a is nuclear removal, f is the fraction of nuclear outflow exported under this model, and b is protein loss. All three remain free. Absolute values of every individual parameter can range over the positive real axis when the others vary. The observations do identify the combinations

    s/(e+d_N)=n, e/d_C=c/n, k/d_P=p/c.

Thus it would be incorrect to say that exact compartment measurements identify no rate comparisons at all. They identify e/d_C, and positivity gives the sharp open bound

    0 < d_C/(e+d_N) < n/c,

whose endpoints are limiting values as f tends to 0 or 1. No positive lower bound on f, no upper bound below 1, and no point estimate of it follows from the equilibrium triple alone.

The common-clock symmetry multiplies all six parameters, including the synthesis flux s, by one lambda>0. It preserves equilibrium and changes only physical speed. It removes only one of the three free dimensions. The export fraction f and relative protein/nuclear clock b/a can still change, so supplying one global time scale is not sufficient.

## Exact synthetic counterexamples

All rows below yield exactly (n,c,p)=(2,3,5). Parameter order is (s,e,d_N,d_C,k,d_P).

| Row | Parameters | Distinction |
| --- | --- | --- |
| A | (8,3,1,2,5,3) | reference |
| B | (8,1,3,2/3,5,3) | same synthesis, same nuclear removal and same protein clock, but export fraction 1/4 rather than 3/4 |
| C | (8,3,1,2,10,6) | only protein production/loss clock changes |
| D | (16,6,2,4,10,6) | common time rescaling by 2 |

A and B are not common rescalings: their shared positive d_P fixes a putative scale factor to 1, yet their e differs. A and C similarly defeat an explanation that only one shared clock is unknown. The example falsifies point identification even with exact amounts and known s and e+d_N.

## Additional information that could separate the rates

This is a statement of information requirements, not a proposal for experiments or interventions.

For exact compartment amounts, one sufficient set of three further scalar constraints is known a=e+d_N, known d_C, and known d_P. Reconstruction is

    e=d_C c/n, d_N=a-d_C c/n, s=an, k=d_P p/c,

with consistency a>d_C c/n. Equality lies on the d_N=0 boundary and is not in the strictly positive model. Equivalently, authenticated synthesis flux s and export flux J=e n, with 0<J<s, plus protein production flux H=k c, determine

    e=J/n, d_N=(s-J)/n, d_C=J/c, k=H/c, d_P=H/p.

Such constraints must supply independent information about nuclear throughput/partition and the separate protein clock. Three independent regular smooth scalar constraints are generically necessary to remove the three-dimensional fiber locally; rank alone is not a proof of global uniqueness for an arbitrary nonlinear constraint system. One known scalar rate, a single generic half-life, further steady-state repeats, or an uncalibrated fluorescence value do not suffice in general.

Under the model, compartment-resolved kinetic information that supplies a and d_C, together with calibrated c/n, separates export from nuclear degradation. A distinct protein turnover or production constraint is still needed for k and d_P. A borrowed rate from another transcript, construct or state is an assumption, not a measurement of this system. Assuming d_N=0 changes the model and fixes f=1; it does not follow from the equilibrium.

## Boundaries and weaker observations

Strictly positive parameters imply strictly positive equilibrium amounts. An exact zero amount is inconsistent with this domain. If zero rates are admitted, new boundary strata arise; for example d_N=0 is regular with e>0, whereas d_P=0 and k,c>0 admit no finite protein equilibrium. If both k and d_P vanish, arbitrary protein levels can be stationary and stability/uniqueness fails. Those boundary strata are not silently included in (1).

Even a perfectly calibrated **total** RNA amount t=n+c and calibrated protein p leave the additional free split 0<n<t, c=t-n. Then (1), with that free n, is a four-dimensional family. Its physical bound k/d_P=p/c>p/t is sharp but depends on calibrated quantities and the model. Deposited standardized assay effects are not those quantities.

A diagnostic pair illustrates a further interpretive error. Let a reference have (n,c,p)=(1,1,1) and a second state (3,1/2,1/2). Both can have exactly k=d_P=d_C=d_N=1, with reference e=1,s=2 and second-state e=1/6,s=7/2. Total RNA rises from 2 to 7/2 while protein falls from 1 to 1/2, without any change in k or protein degradation. Opposite RNA/protein directions alone therefore do not identify translational compensation, even in an idealized calibrated setting.

## Dependency and obligation ledger

| Obligation | Status | Evidence |
| --- | --- | --- |
| All positive parameter sets characterized, both directions | passed | three equilibrium equations and inverse parameterization (1) |
| Fiber dimension and distinction from common scale | passed | positive smooth coordinate bijection; explicit rows A/B/C/D |
| Exact positive/zero boundary scope | passed within stated scope | explicit exclusions and boundary counterexamples above |
| Stability versus actual assay equilibration | conditional for assay | stability proved; assay initial conditions and rates unavailable |
| Additional information sufficiency | passed mathematically | substitution of reconstruction equations |
| Applicability to COMT physical rate inference | inapplicable with current proxy contract | decision_report.md |
| Primary literature comparison | checked and bounded | prior_map.md; no priority claim |
| Universal proof supported by finite tests | not asserted | synthetic checks are arithmetic regression checks only |
