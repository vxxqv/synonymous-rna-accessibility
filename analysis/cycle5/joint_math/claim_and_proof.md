# Joint trait preservation: a bounded coding consequence of Steinitz

Status: author derivation and self-review only. Fresh adversarial review is required before use. No biological outcome has been read. Generic bottleneck optimization, vector balancing, reconfiguration hardness and Laplacian connectivity are established prior work, not new methods here.

## Contract and question

Let G be a finite undirected genotype graph. An edge changes exactly one nucleotide and preserves the specified protein sequence and any separately declared compatibility obligations. All vertices, not just endpoints, must satisfy those obligations. Endpoints s,t belong to G and satisfy every joint phenotype requirement. Each trait has an externally specified positive scale a_k. This is not an inferred biological calibration.

For lower requirements f_k(x) >= tau_k, put

    psi(x) = max_k max(0,(tau_k-f_k(x))/a_k).

For acceptable closed intervals I_k=[L_k,U_k], put

    psi(x) = max_k max(0,(L_k-f_k(x))/a_k,(f_k(x)-U_k)/a_k),

omitting an infinite endpoint. At least one finite bound is assumed per included trait. Let

    R(s,t) = min over finite s-to-t paths P of max_{x in P} psi(x).

Use R=+infinity if s and t have no base-graph path. For s=t, R=psi(s)=0 under the endpoint eligibility condition. Any repeated vertex can be removed without increasing the maximum, so in a finite graph the minimum is attained on a simple path. No phenotype relaxation can reconnect distinct base components. For each epsilon >=0, R<=epsilon exactly when s,t are connected in G[{psi<=epsilon}]. This is a standard vertex-bottleneck path problem.

Trait intervals are important: a quantity with an acceptable optimum or window cannot simply be reoriented into a globally increasing biological benefit. R measures declared assay compromise, not fitness, energy, evolutionary probability or time.

## Main theorem: trait-effect rank plus a local coding excursion

Let G be a Cartesian product of finite connected factor graphs H_1,...,H_m, with actual edges changing one factor along an edge of that factor. Let the scaled trait vector have an additive representation

    F(x)= (f_1(x)/a_1,...,f_q(x)/a_q)
        = c + sum_i phi_i(x_i).

Let K=product_k [L_k/a_k,U_k/a_k], with lower/upper-only intervals permitted, and assume F(s),F(t) belong to K. Let D={i:s_i!=t_i}; put n=|D|. For i in D define

    v_i = phi_i(t_i)-phi_i(s_i),
    M = max_i ||v_i||_infinity,
    r = dim span{v_i:i in D}.

For each changed factor, define its endpoint box B_i by the coordinatewise minima and maxima of phi_i(s_i),phi_i(t_i). Define

    eta_i = min over H_i paths Q from s_i to t_i
            max_{z in Q} dist_infinity(phi_i(z),B_i),
    eta = max_i eta_i.

For n=0 set r=M=eta=0. Then there exists a nucleotide-valid factor-block route with

    R(s,t) <= r M + eta.                                      (C5.1)

Its changed factors are processed once each, although a path inside a factor need not be a single direct nucleotide change. The route need not improve every trait monotonically. If all endpoint scaled values have uniform distance at least sigma>=0 from every finite face of K, the sharper margin statement is

    R(s,t) <= max(0,r M+eta-sigma).                            (C5.2)

The factor r is a safe classical constant, not an asserted best Steinitz constant. The rank term is attained in a rank-one binary example below; the coefficient of eta cannot be reduced in general. Neither establishes optimality of r in higher rank. This theorem is a coding/product-graph corollary of classical vector balancing, not a claimed new Steinitz theorem or an exact algorithm for R.

### Proof

The case n=0 is immediate. If M=0, every block boundary has phenotype F(s)=F(t), irrespective of block order. The local-box argument below gives the bound eta, so suppose M>0 and r>=1.

Work in the r-dimensional real subspace V spanned by the v_i, with the restriction of the ambient infinity norm. Let d=sum_i v_i. Add n labelled dummy vectors, each equal to -d/n. Every dummy has norm at most M, because ||d||<=nM, and the augmented list sums to zero. The classical Steinitz lemma in an arbitrary r-dimensional norm gives a permutation whose every partial sum has norm at most rM. The external leaf is exactly this zero-sum statement, as given in Barany2024 Theorem1.2 and attributed there to Grinberg-Sevastyanov1980.

Delete the dummy entries, keeping the induced order of the original factors. Immediately after the k-th original vector in the augmented list, let j be the number of dummies already seen. Its partial sum is

    sum_{h<=k} v_pi(h) - (j/n)d.

It has norm at most rM and 0<=j/n<=1. Thus every block-boundary phenotype lies within rM of the endpoint segment [F(s),F(t)]. This argument uses one common scalar j/n for the entire vector, so mixed signs in d cause no problem. The segment is contained in K, since K is convex. Consequently every boundary belongs to K+[-rM,rM]^q. The zero boundary is F(s); the last boundary is exactly F(t).

Choose each factor path to attain eta_i, which is possible in a finite connected graph after removing repeated vertices. During block i, every coordinate of the phenotype lies between the corresponding before-block and after-block coordinates, up to error eta_i: this is exactly the definition of distance to the local endpoint box B_i. Both boundary coordinates lie in the corresponding expanded acceptable interval. Hence all block intermediates lie in K+[-(rM+eta),rM+eta]^q. This proves C5.1. This lifting step needs an axis-aligned endpoint box, not proximity to the much smaller endpoint line segment.

If the two endpoints lie sigma inside every finite face, their entire segment lies in the shrunken box with the same margin. Repeating the argument gives C5.2. Negative resulting violation bounds are clipped at zero. There is no claim of a calibrated value of sigma.

### Application to exact standard-code synonymous graphs

For one fixed, nonoverlapping coding frame with standard sense codons, the unrestricted synonymous nucleotide graph is the Cartesian product of the single-codon synonymous graphs. Restrict to the component containing s: every factor has at most six codons. Serine's UCN and AGY groups are different components. If t is in another such component, R=+infinity. Other amino-acid factors are connected. Fixed bases simply induce factor subgraphs, and one must first take their actual connected components.

An edge within a codon factor is a one-nucleotide substitution. Leucine and arginine factors are not cliques; a direct replacement of their endpoint codons can require intermediary synonymous codons. Eta handles that distinction exactly. It is a bottleneck query on a graph of at most six vertices, with rational vertex costs when input effects/scales are rational. A product of codon factors is not generally a median graph: fourfold wobble factors contain triangles. No median assumption is used.

RNA pairing constraints between codons, overlapping translated frames, motifs spanning factor boundaries, or other coupled obligations can destroy this product. C5.1 is not asserted for them. A factorization into larger genuinely independent blocks is permitted, but factor sizes and local bottleneck costs must then be accounted for explicitly.

## Quantitative interpretation and boundaries

The result controls a model-based maximum temporary compromise by intrinsic effect rank and the largest scaled codon change, rather than by the number of changed codons. Adding redundant measured traits does not automatically increase r. It can still change M, eta, margins and actual acceptability. Additive representation, scales and coverage of local intermediates are substantive assumptions. An empirical obstruction above this bound would refute at least one of those assumptions or the trait estimates; it would not by itself identify an RNA mechanism.

**Deterministic approximation corollary.** Suppose true scaled traits F_true need not be additive, but ||F_true(x)-F_add(x)||_infinity <= epsilon at every relevant base-graph vertex, where F_add has the representation above. If the true endpoints belong to K, then

    R_true(s,t) <= r M + eta + 2 epsilon.                     (C5.3)

Indeed the model endpoints belong to the expanded box K^epsilon. Apply C5.1 to that box, then apply the uniform error bound once more along the resulting route. Expansion of a rectangular box adds the radii. This requires a global deterministic bound over all relevant states, not training residuals, a held-out mean error, correlation, or an uncalibrated confidence interval. The coefficient2 cannot generally be decreased: on the Leu example let F_add be the constant -epsilon, true endpoints0 and every other true codon -2epsilon. Then r=M=eta=0, the uniform error is epsilon, but R_true=2epsilon. No calibrated epsilon is assumed or recovered here.

For lower bounds only, a simpler formal corollary is useful. Suppose every changed block's vector v_i is either componentwise nonnegative or componentwise nonpositive, and its chosen local path is componentwise no worse than the lower of its two local endpoints. Process the nonnegative blocks first and the nonpositive blocks second. Every intermediate trait is at least the lower of its two global endpoint values. Thus R=0. This is a shared-order sufficient condition, not necessary. It fails for interval requirements because an early beneficial overshoot can cross the upper bound. It also fails if local coding gates violate the stated local-path condition.

## Counterexamples that must remain attached to the claim

1. **Conflict with no epistasis.** On Q_2, f=(x_1-x_2,-x_1+x_2), thresholds(0,0), s=00,t=11. Each trait separately permits a path, but jointly both intermediates violate one trait by1. R=1, rank=1,M=1,eta=0. This is only a minimal illustration and rank-one sharpness check, not a novelty claim. The two forbidden intermediates form an irredundant interior separator, and neither trait alone creates that obstruction.

2. **Local coding excursion is essential, even at rank zero.** One Leu codon, s=UUA,t=CUC. Set a scalar contribution0 at both endpoints and -3 at the other four Leu codons, with requirement f>=0. They are not adjacent. Every path has an interior codon, all such codons have score -3. Therefore R=eta=3, while every endpoint effect is zero, so r=M=0. Any claim R<=rM without a codon-path correction is false.

3. **A block-order search need not find the global optimum.** On Q_3 set f=(x_1-x_2+x_3,-x_1+x_2+x_3), thresholds0, s=000,t=110. Every ordering of the two endpoint changes has compromise1. The simple route 000,001,101,111,110 has compromise0 by temporarily using a background coordinate. Thus global R is0, and restricting search to endpoint-difference permutations can give a false obstruction. The main bound provides an upper bound, not an assertion that its block witness is optimal.

4. **The product assumption is essential.** In the independent-set flip graph of K_{2,2}, let f(I)=|I| and require f>=2. Take the two full opposite parts as endpoints. Any route must pass through the empty set, so R=2. Treating its four binary variables as uncoupled factors would give r=1,M=1,eta=0 and the false bound1. This constrained graph occurs as a fixed-orientation simultaneous-pairing logical graph in the existing code-specific theory. It is not a claim that one physical RNA structure contains all those pairs at once.

5. **Closed versus open requirements.** On the Leu example give endpoints score1 and all other codons score0. The weak requirement f>=0 connects them and gives R=0. The strict requirement f>0 disconnects them. Positive-part violation loses strict-boundary information. For strict thresholds use the signed violation h=max_k((tau_k-f_k)/a_k), or its interval analogue. With T=min_P max_v h(v), strict relaxed connectivity h<epsilon holds exactly for epsilon>T, not at equality. Finite T is attained as a path value, but the open-filtration set of sufficient epsilon has an unattained infimum T.

6. **Baseline disconnection.** Ser UCU and AGU encode the same amino acid but lie in different synonymous components. Every R is+infinity irrespective of phenotype values.

7. **Interval alignment is insufficient.** Q_2, scalar f=x_1+x_2, interval[1,1], endpoints10 and01. Both codon changes have scalar effects, but any route leaves the equality slice. R=1 with unit scale. Lower-bound shared-order connectivity must not be extended to arbitrary acceptable intervals.

## Complexity and bit conventions

For an explicitly listed finite graph, evaluating q rational traits costs O(q|V|) arithmetic operations; a minimax Dijkstra implementation or critical-value activation solves the remaining standard bottleneck query in polynomial time in graph size. Rational comparisons and bit lengths remain part of bit complexity. This graph can be exponential in coding length. Paths can have exponential length in a succinctly specified instance, even though a loop-free path has at most |V|-1 edges.

The rank and M are computable by exact rational linear algebra from the endpoint-effect matrix. Every eta_i uses at most six codon vertices; costs are coordinatewise rational box distances. The written existential proof imports classical Steinitz and does not implement or claim a new polynomial-time ordering algorithm. Exhaustively trying n! block orders is a finite witness method, not an efficient general algorithm.

Cardinal et al.2018 Theorem4.7 already makes exact joint reachability PSPACE-complete with a constant number of affine requirements on Q_n. Their displayed construction uses70 dimensions with lower/upper constraints. Repeating the standard two-codon Lys choice AAA/AAG realizes Q_n exactly while preserving one peptide. Thus the prior hardness transfers directly to this biological graph class, even with constant trait rank bounded by70. Do not claim hardness for two traits from that result. Fixed trait count, fixed rank, additive scores or a sparse interaction description alone is not a tractability theorem. The relaxation bound does not conflict with exact hardness.

## Self-audit and dependency graph

Definition of R -> finite simple-path reduction and bottleneck equivalence: internally proved, passed.

Product graph plus additive traits -> endpoint effect decomposition: assumption, passed as typed; biological validity unestablished.

Classical zero-sum Steinitz in the intrinsic subspace -> augmented partial-sum bound: external result verified in Barany2024 Theorem1.2; norm/rank translation passed.

Collinear dummy removal -> proximity to the endpoint segment: internally proved, passed including mixed-sign d.

Endpoint box plus exact factor path -> full nucleotide lift: internally proved, passed; eta cannot be omitted.

Computational fixtures -> implementation checks only in their declared range, not proof of C5.1.

Primary self-review verdict: **conditional on the named classical Steinitz input**, with a complete internal reduction and lifting proof. No independent reviewer has yet audited this new package. Closest-prior classification: **formal corollary not stated** for the rank-bound/codon-lifting formulation, and **known after translation of notation** for bottleneck evaluation and the fixed-trait hardness. Global novelty is not established.
