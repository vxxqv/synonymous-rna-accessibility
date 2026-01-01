# Independent review of the equilibrium-observation and processing claims

**Verdict: the global positive-domain classification and the conditional processing-identifiability argument pass.** All 1,024 observation subsets and all 51 minimal sufficient triples have inspectable global evidence. No positive-domain counterexample to the stated classification was found. One conditioning sentence requires more precise terminology; no rate classification changes.

This is a separately authored review with inherited project context, not a blinded external review. The reviewer inspected the claims and source code but did not import or execute the author's mathematical checker. Only synthetic expressions, abstract mixtures, and an isolated source parser function on symbolic slice predicates were evaluated. No source phenotype rows, raw reads, new nucleotide sequences or fitted outcomes were used.

## G: independent global argument

Use different coordinates from the author's checker: total nuclear exit a=e+d_N, export fraction f=e/a and protein loss z=d_P. Their domain is a,z>0 and 0<f<1. For fixed calibrated n,c,p>0, the full equilibrium fiber is

```
s=n*a, e=a*f, d_N=a*(1-f), d_C=n*a*f/c, k=p*z/c, d_P=z.
```

The three balance equations hold identically, and the coordinate change is bijective. The observation functions reduce, after known positive calibration factors, to a, a*f, a*(1-f), f and z.

Two distinct nuclear quantities have the following global inverses. The last column defines the positive image, not a numerical filtering rule.

| Known normalized quantities | a | f | Admissibility |
|---|---|---|---|
| A,E | A | E/A | A>E>0 |
| A,D | A | 1-D/A | A>D>0 |
| A,F | A | F | A>0, 0<F<1 |
| E,D | E+D | E/(E+D) | E,D>0 |
| E,F | E/F | F | E>0, 0<F<1 |
| D,F | D/(1-F) | F | D>0, 0<F<1 |

One normalized protein observation gives z. There is no second positive branch, sign choice or excluded denominator in these domains. A larger observation set must satisfy all redundancy and cross-class consistency equations. An inconsistent vector has no solution; it is not an identified admissible state.

The independent symbolic checker considered all 120 triples. A Jacobian determinant was used only to nominate candidates, not to establish global sufficiency. Each of the 51 candidate triples was then solved symbolically, and both compositions of its inverse and observation map were verified as rational identities. The image conditions above complete the global proof. This is stronger than the author's finite five-anchor checks and does not depend on them.

For every insufficient subset, at least one of the following globally positive, nontrivial maps leaves all its selected observation functions unchanged:

| Map in (a,f,z) coordinates | Observations it can preserve |
|---|---|
| (a,f,2z) | All nuclear quantities; no protein clock |
| (2a,f,z) | f and the protein quantities |
| (a,f/2,z) | a and the protein quantities |
| (2a,f/2,z) | a*f and the protein quantities |
| (2a,(1+f)/2,z) | a*(1-f) and the protein quantities |

Every map stays inside the strict domain at every starting point. The checker computes invariant observation masks directly from these symbolic functions and finds an appropriate map for every one of the 219 insufficient subsets. It does not assume the author's class-based Boolean classifier. Every other subset contains a verified inverse triple. This exhausts all 1,024 subsets and confirms 805 sufficient subsets and 51 inclusion-minimal triples.

The reasoning is global because it supplies inverses on the full positive image or a distinct positive alternative at every point, not merely generic rank. The fixed candidate list is essential. The count cannot be promoted to a statement about arbitrary assays, trajectories, nonlinear measurements or a broader kinetic model.

## Occupancy and residence time

Under the stipulated ideal relation O=p*z*tau, the map z -> 2z, tau -> tau/2 preserves occupancy and every nuclear observation. The four nuclear ambiguity maps above preserve occupancy with tau unchanged. Thus an unknown residence time prevents O from replacing a protein-clock observation or repairing an insufficient nuclear set.

If the original rates are identified and O is known, tau=O/(p*z) is uniquely positive. Both inverse compositions were checked symbolically. Without O, changing tau alone leaves every original candidate unchanged. Consequently the augmented rate-plus-residence vector is identified exactly when O and an originally sufficient set are observed. The 51 minimal triples become 51 minimal four-observation sets after appending O. The optional complete enumeration over the 2,048 original-plus-O subsets gives 805 augmented sufficient sets; identifying rates alone remains possible in 1,610 sets because O can be present or absent without changing that criterion.

This conclusion depends on the asserted O=H*tau relation. It does not establish that a heavy-polysome bin measures O, that tau is known, that occupied ribosomes finish proteins, or that initiation equals completion. Export flux J and completed-protein flux H remain distinct.

## Boundary and calibration challenges

The strict-domain qualification is necessary. If d_N=0 and f=1 are admitted, the supposedly sufficient positive-domain triple {d_N,f,d_P} takes values (0,1,1) for every positive e with d_P=1. The positive-domain inverse D*F/(1-F) is undefined there. This does not falsify the theorem, which excludes that boundary, but prevents silent extension to it. A model with d_N=0 known in advance has a different sufficient-observation problem, as the author states.

Calibrated amounts are also necessary. An explicit outside-contract example retains observed n_tilde=c=p=1 and observed scalars s=4,e=1,d_P=1:

- Physical n=1 and nuclear gain 1 give d_N=3,d_C=1,k=1.
- Physical n=2 and nuclear gain 1/2 give d_N=1,d_C=2,k=1.

Both satisfy all equilibrium balances and have the same displayed amounts and candidate observations, but different rates. The triple {s,e,d_P} is sufficient only when the physical n is authenticated, not when an unknown gain turns a relative readout into an apparent amount. The example is synthetic, not a claim about actual source measurements.

The model also requires compatible steady state, a common population and time scale, and consistent effective dilution/loss. Snapshots from different compartments or selected species cannot be substituted merely because their labels resemble n,c,p.

## Required conditioning wording

For e=D*f/(1-f), the unambiguous derivative is

```
d(log e)/df = 1/[f(1-f)].
```

This is relative output sensitivity per absolute change in f. The relative-input elasticity is instead

```
d(log e)/d(log f) = 1/(1-f).
```

Both were verified symbolically. Replace “relative fraction sensitivity” with the first explicit derivative, as root proposed. Near-cancellation in A-E and the fraction-boundary conditioning remarks otherwise hold. Structural uniqueness must not be presented as precise practical estimation.

## C: known-matrix statement and positive mixture counterexample

For a known finite matrix A on a full-dimensional open positive domain, a linear functional l*x is globally determined by A*x if and only if l lies in the row space of A. Sufficiency is immediate from l=b*A. If l is outside the row space, some kernel vector w satisfies l*w != 0. At any interior x, sufficiently small perturbations x +/- epsilon*w are positive and have the same A*x but different l*x. An observation on the nonnegative boundary can identify an otherwise non-estimable contrast accidentally; that does not contradict this global open-domain statement.

For the author's abstract aggregation matrix

```
A = [[1,0,0,0],
     [0,1,1,0]],
```

the rank is two and its kernel has the stated two dimensions. The coding contrast (1,1,0,0) is outside the row space even after adding total amount. An independent positive example is u=(3,1,4,2), v=(3,3,2,2). Both have observed counts (3,5) and total 10, while the proposed coding contrast is 4 versus 6. The author's separate example, totals 7 and coding values 5 versus 3, also checks exactly.

All 16 binary contrasts were independently classified. Without total, a contrast is identified precisely when its two hidden-junction coefficients agree and its unretained-species coefficient is zero. With total, equality of the two hidden-junction coefficients remains necessary and sufficient. Exact positive kernel perturbations certify every nonidentified case.

This is conditional linear estimability, not an empirical rank calculation for Dao. The actual transcript-to-category recovery and coding-product maps are unauthenticated. If the actual coding label were constant within each reported category, that particular hidden-path ambiguity would disappear. Conversely, unknown recovery can add ambiguity. The author appropriately treats the example as a possible unresolved refinement rather than asserting actual hidden isoforms.

## Source-only parser and median checks

Source `categorize_read`, lines 30-47, accumulates recognized deletion lengths while overwriting donor/acceptor variables at each recognized deletion; lines 65-71 report the last such junction. Isolated execution of this original function on two symbolic CIGAR/predicate fixtures confirms identical final category output despite an additional earlier deletion. No nucleotide reference was synthesized and no real read or path occurrence was claimed.

The inspected processing code independently supports the stated operator:

- `process_raw_count`, lines 4-23, retains full-length/spliced rows, rejects retained reporter totals <=10, then divides by retained reporter totals.
- `filter_barcode_count`, lines 25-31, counts full-length rows per RE, not distinct barcode identities. Its threshold 5 retains counts strictly greater than 5, so three barcodes repeated across two preparations can pass as six rows.
- `compute_median_splice_fraction`, lines 33-40, calls a pivot without an explicit aggregation function, then takes coordinatewise medians. Missing category cells are filled after aggregation.

The author's six-barcode example yields medians (0.1,0.1,0.1), total 0.3. A different exact synthetic example using repeated normalized vectors (0.6,0.3,0.1), (0.1,0.6,0.3), (0.3,0.1,0.6) yields (0.3,0.3,0.3), total 0.9. Both use six barcodes, 100 counts each and positive full-length counts. The present-only pivot example gives full-length 0.95 plus spliced 0.1, rather than the explicitly zero-completed spliced mean 0.05. These computations verify nonlinearity and missing-cell semantics, not actual source artifact rates or invalidity of the source's intended statistic.

No full historical pandas pipeline was replayed. This review checks the source-defined operator and exact invented arithmetic; it does not authenticate a historical runtime or hidden source-data distribution.

## Supported increment and application limits

The supported mathematical increment is a complete elementary global classification for this specified equilibrium fiber and ten candidate scalar observations, together with explicit missing-information examples. Equilibrium elimination, global identifiability, linear estimability and observation aggregation are established ideas. The review supplies no broader priority claim.

Neither COMT's low-bin enrichment/load categories nor cPRIME's local RNA/DNA representation is authenticated as the calibrated n,c,p or kinetic observations required here. Neither supplies a known residence time or a complete transcript-to-product observation matrix. These theorems can delimit what additional records would be required; they cannot turn those current proxies into identified translation, degradation, export or protein-product mechanisms.

## Reproducibility and immutable review point

The reviewed rate claim is SHA-256 `72d2b9d519ac375b55a70ca573501b53bcf440dc6fb9bfbd8df1556aa94af3dc`; the processing claim is `f2a387383c8ed130b1b4f6a846b07c9364bc8eb182668ea578c3165339145ee2`. The source hashes and every symbolic inverse, subset proof type, counterexample family and processing check are saved in `exact_v2/`.

`verify_symbolic.py` uses CPython 3.12.14 and the already installed SymPy 1.14.0. Its 2,497 checks include rational-function identities rather than a new set of numerical rate anchors. It refuses overwrite and optimized execution and rechecks source hashes. Run with explicit `--source-dir`, `--dao-gate` and a fresh `--out` directory. No network or Git is needed.

The first review draft is retained in `history_v1/verify_symbolic.py` and `exact_v1/`. Before handoff, its occupancy inverse check was strengthened into two explicit compositions plus the missing-O ambiguity. `exact_v2/` contains the authoritative receipt. All scientific evidence outputs are byte-identical between these review runs; no author input or claim was changed.
