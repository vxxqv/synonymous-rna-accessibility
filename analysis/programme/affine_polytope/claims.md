# Multidimensional correlated uncertainty: reviewed formulation

Status: separately authored, non-blinded proof and exact computational review completed. This second revision makes the disclosed wording and complexity corrections below. These are abstract finite graph statements. No biological uncertainty region is calibrated here. Baseline: repository e9846dca9f74b2c5c34e984eed970e996b9df2bf and the one parameter affine result in the current manuscript.

## Definitions

Let G be a finite undirected graph with n vertices and m edges. Distinguished vertices s,t are connected; paths include endpoints and may be simple. For a field f, C(f) = max_P min_(v in P) f_v, e(f) = min(f_s,f_t), and B(f) = e(f) - C(f). A singleton path is allowed when s=t. An empty compatibility path set is reported separately with extended capacity minus infinity and barrier plus infinity.

Let Q = {theta in R^d: A theta <= b} be a nonempty bounded rational polytope, with k inequalities. Scores f_v(theta) = a_v + q_v dot theta have rational coefficients. Q need not have full dimension. Fixed, guaranteed and possible capacities have respectively the quantifiers exists P for all theta, for all theta exists P, and exists theta exists P.

## U1. Exact fixed dimension identification

Form a list H of the k boundary hyperplanes of Q and every nonconstant score equality f_i = f_j. Ignore identities and inconsistent constant equalities. For each d element subset of H with linearly independent normals, solve their intersection and retain it if it belongs to Q. Remove duplicate points. For d=0 use the single empty coordinate vector. Call the retained set T.

The minimum and maximum values of C and B on Q are attained at points in T. This does not assert that every optimizing point lies in T. Consequently exact capacity and barrier ranges require at most binomial(k + n(n-1)/2,d) rational candidate solves and widest path evaluations. Coordinate minima over Q give the fixed route capacity through one additional widest path calculation. The guaranteed and possible capacities are the minimum and maximum of C over T. No parameter corner shortcut is valid in general.

Proof. The score equality arrangement partitions Q into finitely many closed polyhedral cells with fixed weak score order. For each cell choose any point in its relative interior. Every path bottleneck is a fixed affine score on that relative interior, modulo identically tied scores. An optimal path there stays optimal throughout the relative interior. Thus C is affine on the cell; continuity extends the formula to its boundary. The smaller endpoint score is affine on the cell, so B is affine there as well. A nonempty bounded cell is the convex hull of its vertices. Each vertex has d linearly independent active constraint normals in the ambient space, even for a lower dimensional cell. All these boundaries belong to H. Thus it is in T, and both extrema of each affine cell function occur in T. Coordinate minima also occur at Q vertices, which are included in T. The fixed route identity is max_P min_v min_theta f_v = C(ell).

Arithmetic bound for the elementary implementation: with H = k + n(n-1)/2, constructing, normalizing and sorting the possible boundaries requires O((H+1)(d+1) log(H+2)) arithmetic/comparison work, including score-pair preprocessing when d=0. After this preprocessing there are at most binomial(H,d) candidates, each costing O(d^3 + kd + nd + (n+m) log(n+1)) arithmetic/comparison operations. Dependent intersections are discarded. Exact rational bit complexity is additional but polynomial in input bit length for fixed d, by rational linear system bounds and polynomially many arithmetic operations on bounded determinant sizes. This is polynomial for each fixed d, not a fixed parameter tractable bound and not polynomial uniformly in d. The implementation initially accepts a bounding box plus optional rational inequalities so boundedness is explicit.

The arrangement idea and lattice polynomial structure are established tools. This result is an explicit accessibility application and extension of the existing one parameter statement, not a priority claim for general piecewise linear optimization.

## U2. A variable dimension boundary

With d part of the input, deciding whether the possible capacity is at least 1 is NP complete, even on a chain of clause layers of width at most three, with theta in [0,1]^d and scores only 1, theta_i and 1-theta_i.

Reduction. From a 3CNF formula, create one layer per clause and one vertex per literal occurrence. Join all vertices in consecutive layers, add s joined to the first layer and t to the last, and give s,t score 1. A positive occurrence of variable i has score theta_i; a negative occurrence has score 1-theta_i. Any s to t path must encounter every layer. At threshold 1 it therefore selects at least one true literal per clause, with selections consistent because each parameter has one value. Parameters that appear in selected vertices are forced to 0 or 1; extend the remaining coordinates to a Boolean assignment arbitrarily. Conversely, a satisfying assignment supplies one retained vertex per layer and hence a path. The graph is undirected, which does not affect the layer separation argument. A satisfying path and feasible rational parameter are a polynomial certificate; for general rational affine input, feasibility of its inequalities has a polynomial bit rational witness.

Deciding whether the guaranteed capacity is at least 1/2 is coNP complete, even when the graph consists of internally disjoint clause paths of at most three internal vertices sharing only s,t, and with the same parameter box and score alphabet.

Reduction. For each clause create a separate path from s to t with one internal vertex per literal. Give a positive literal vertex score 1-theta_i and a negative literal vertex score theta_i. A clause path is retained at threshold 1/2 if all of that clause's literals can be false under the parameter values. If the formula is satisfiable, its Boolean satisfying assignment makes each clause path contain a score zero, so no path remains. If it is unsatisfiable, round each theta_i to a Boolean value, declaring true when theta_i >= 1/2. Some clause is false under this assignment. Every positive literal in that clause has theta_i < 1/2, hence complement score > 1/2; every negative literal has theta_i >= 1/2, hence score >= 1/2. Its entire path remains. Thus every parameter scenario connects exactly when the formula is unsatisfiable. Ties satisfy the closed threshold convention. For membership, a failure certificate consists of a vertex separator K and rational theta such that f_v(theta) < threshold for all v in K. Strict feasibility can be represented by a positive common slack variable; rational linear programming gives a polynomial bit positive witness whenever one exists. A separator may include s or t. Graph separation and rational strict inequalities can be checked in polynomial time.

These reductions address a shared parameter representation of scores. An independent interval for every vertex remains easy. No hardness claim is made for actual RNA measurement models or fixed latent dimension. Exact prior overlap is unresolved until a dedicated search and review are complete.

## U3. Structural classification of route adjustment

A minimal path set is an inclusion minimal vertex set whose induced subgraph connects s to t. Let M be the finite family of such sets. The following are equivalent:

1. M contains exactly one set.
2. For every nonempty compact uncertainty set of vertex scores, fixed and guaranteed capacities coincide.
3. The same equality holds for every one parameter affine score family on [-1,1].

Proof. If M={P}, every connecting vertex set contains P, and C(f)=min_(v in P) f_v. The two minima commute for any uncertainty set. The second statement implies the third. For the converse choose distinct minimal sets P,Q. Their intersection contains no connecting path, since otherwise P or Q would fail minimality. Put score 1 on P intersection Q, theta on P minus Q, -theta on Q minus P, and -1 outside their union. For theta>=0, P supplies capacity at least zero; for theta<=0, Q does. At theta=0, their intersection is the only score 1 set and is disconnected, so the guaranteed capacity is zero. The coordinatewise lower field has score -1 off P intersection Q and 1 on it, hence the fixed capacity is -1. This contradicts the third statement.

A direct s,t edge yields the unique minimal vertex set {s,t}, even when other simple paths exist. Thus uniqueness of a simple path is sufficient but not necessary. The result uses standard minimal path sets and lattice polynomial representation. Whether this exact quantifier characterization is explicitly established elsewhere remains to be checked.

## Completed review and remaining limits

The separately authored review challenged lower dimensional polytopes, ties, strict versus closed thresholds, certificate bit lengths, undirected backtracking and minimal vertex path sets with chords. An exact Fraction/Fourier-Motzkin path and separator implementation agreed on 81 original case statuses and 395 rational quantities. The review was not blinded and inherited project context. Its proof arguments, receipt and source hashes are retained in ../independent_review/.

The original case named two_dimensional_corner_failure had C=|theta1+theta2| and already attained its extrema at square corners. Its original input and answer are preserved; this revision calls it rank_one_sum_two_parameters. A separately disclosed six-vertex case with four parallel intermediate scores theta1, -theta1, theta2 and -theta2 has C=max(|theta1|,|theta2|) and endpoint scores one. All four square corners have C=1 and B=0, while the origin has C=0 and B=1. The reviewer checked this author-suggested challenge separately before this revision.

Marichal's weighted lattice-polynomial normal forms and score-order representation are established prior. U1 is an accessibility application of that machinery; U3 is a short structural consequence. U2 has explicit reductions and membership proofs, but exact priority against correlated-network optimization remains unresolved. No generic algorithm or first hardness result is claimed. Finite checks never prove U1-U3. No biological validation or statistical coverage follows from these deterministic results.
