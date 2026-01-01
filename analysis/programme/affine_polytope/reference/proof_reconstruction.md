# Independent reconstruction before consulting recorded outcomes

This is a separately authored review with inherited project context. It is not blinded, and no new isolated reviewer was available. The reviewer did not author U1-U3. This note was written before reading the author's recorded result receipt. The author's solver was read for implementation review, but will not be imported by the independent computation.

## U1

For a fixed total weak order of the vertex scores, the bottleneck value of every path is one coordinate of that order. Taking the maximum over the finite set of simple paths therefore selects a fixed coordinate on the whole order region. Tied coordinates agree there, so selecting either representative has no effect. The endpoint minimum is also a fixed coordinate. Pulling an order region back along the affine score map and intersecting it with Q gives a bounded polyhedron. Capacity and barrier are affine on it.

At a vertex of a bounded polyhedron, the active ambient constraint normals span the ambient space. Otherwise a nonzero orthogonal direction permits a sufficiently small movement in both directions while preserving inactive strict inequalities and all active equalities, contradicting extremality. This argument also covers a line, a point, or any other rank deficient polytope represented by inequalities. Consequently a vertex is the intersection of d independent boundary or score equality hyperplanes. The extrema are attained among the enumerated points. The statement concerns attainment of extremal values, not a claim that every minimizer belongs to the finite set.

The algorithm is polynomial for fixed d, with exponent depending on d. Rational determinant bounds give polynomial bit lengths for candidates. Generic piecewise affine optimization and lattice polynomial score order formulas are established, so this is an application, not a new generic method.

There is a minor complexity accounting issue in the implementation at d=0: it still constructs all n choose 2 pair equalities. The stated per candidate bound does not count this O(n^2) preparation. A preprocessing term or an early d=0 return resolves that implementation accounting; no extremal value is affected.

## U2

For possible connectivity, a path in the undirected clause chain must cross every layer. At threshold 1 every selected literal fixes its parameter to its satisfying Boolean value. No path can contain selected positive and negative occurrences of the same variable. Backtracking adds restrictions and cannot skip a layer. A Boolean satisfying assignment conversely leaves at least one vertex in every layer and hence connects. The construction is polynomial even with complete bipartite joins of consecutive width three layers.

For guaranteed connectivity, the internally disjoint clause paths mean that any simple endpoint path uses exactly one clause. At a Boolean satisfying assignment, every clause path has an internal zero score and is cut. If a parameter vector disconnected an unsatisfiable formula, rounding theta_i >= 1/2 to true would nevertheless give a false clause. Its positive literals have complement scores strictly above 1/2; its negative literals have scores at least 1/2. That clause path would remain, a contradiction. Values exactly 1/2 cause no failure of the closed threshold convention.

Membership requires proof, not just hardness. For possible connectivity, guess a simple path and a rational point satisfying its closed threshold inequalities and Q. A feasible rational LP has a polynomial bit witness. For the complement of guaranteed connectivity, guess a vertex separator K and a feasible point with f_v < threshold for every v in K. Introduce 0 <= epsilon <= 1 and f_v <= threshold - epsilon. If a strict witness exists, this bounded rational LP has a positive optimum with polynomial bit size and a rational optimizer. Endpoints may be in K. Both membership arguments are polynomial in the explicit graph, rational constraint description, and bit length. They do not establish arbitrary dimension tractability.

## U3

Every connecting vertex set contains an inclusion minimal connecting subset. Thus C(f) is the maximum, over minimal path sets M, of their minimum score. A unique minimal set makes fixed and guaranteed capacity equal by commuting two minima.

If P and Q are distinct minimal sets, neither contains the other. Their intersection cannot connect, since a connecting intersection would be a proper connecting subset of at least one. Scores 1 on the intersection, theta on P minus Q, -theta on Q minus P, and -1 elsewhere give capacity at least 0 for every theta in [-1,1]. At theta=0 only the intersection has positive score, so capacity is exactly 0. The coordinatewise lower scores equal -1 outside the intersection, so fixed capacity is -1. This establishes necessity even with one affine parameter. A direct endpoint edge produces unique minimal set {s,t} despite other longer paths; chords must be handled through minimal vertex sets rather than simple path counts. For s=t the unique minimal set is {s}, and the same equality holds.

## Independent finite representation

The oracle will enumerate minimal connecting vertex sets P and minimal separators K directly from vertex subsets. It will solve rational linear programs by Fourier-Motzkin elimination, without score order cells or importing the author's solver. With E={s,t}:

* max C = max_P max_theta min_(v in P) f_v, a hypograph LP per P.
* min C = min_K min_theta max_(v in K) f_v, an epigraph LP per K.
* min B = min_P min_(j in E) min_theta max_(v in P) (f_j-f_v).
* max B = max_K max_theta min_(j in E,v in K) (f_j-f_v).
* fixed C = max_P min_(v in P) min_theta f_v.

The separator identity follows directly from retained connectivity: C is at most the largest score in every separator; deleting all vertices with score at most C is a separator because the strict superlevel graph cannot connect. The minimum such separator maximum therefore equals C. This argument includes endpoint separators and tied scores.

Finite agreement will corroborate implementation, not prove the universal claims.
