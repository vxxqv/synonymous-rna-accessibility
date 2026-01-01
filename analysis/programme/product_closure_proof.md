# Additive product threshold queries

Status: direct proof with a separate executable implementation, bounded finite checks and an independent derivation and challenge audit. Structural composition of additive merge trees is established prior work. The novelty of the explicit query algorithm below remains unresolved.

## Setting

Let G_i=(V_i,E_i), for 1 <= i <= d, be nonempty finite undirected graphs. Let f_i be a real score on V_i, and let G be their Cartesian product: one move changes exactly one coordinate along an edge of its factor. Write F(z)=sum_i f_i(z_i). Fix a threshold t and a retained starting vertex x with F(x)>=t. Let C_t(x) be its connected component in the induced graph on F>=t. Set N=sum_i |V_i| and M=sum_i |E_i|. The product itself may have product_i |V_i| vertices.

Define A_i initially as {x_i}, and p_i=max_{v in A_i} f_i(v). A coordinate expansion replaces A_i by the connected component of x_i in

    G_i[ f_i >= t - sum_{j != i} p_j ]

and updates p_i. Process coordinates until no A_i changes. Every threshold comparison is inclusive. Choose the lowest vertex index among tied maxima when a canonical maximizer is needed.

## Proposition 1: exact projection closure and connectivity signature

Every expansion enlarges its coordinate set. The algorithm terminates and its final sets satisfy

    A_i = projection_i C_t(x).

There is a reachable vertex z with z_i a canonical maximizer of f_i on A_i for every i simultaneously. Consequently F(z)=sum_i p_i is the maximum score in C_t(x).

For retained endpoints x and y, they are connected at threshold t if and only if their final coordinate sets agree for every factor. Comparing only whether y_i belongs to the closure from x is insufficient.

### Proof

Maintain three invariants after every expansion: each A_i is connected and contains x_i; every value in A_i is at least t minus the sum of the current other-coordinate peaks; and the joint peak anchor z is reachable from x at threshold t. Initially each set is a singleton, the second assertion follows from F(x)>=t, and z=x.

Suppose coordinate i is expanded. Peaks in other coordinates have never decreased, so the new factor threshold cannot exceed the threshold at which the old A_i was obtained. Thus old A_i lies in the new component. Its previous maximizer also lies in that component. Holding every other coordinate at its current peak anchor, any path inside this factor component stays above the product threshold. The anchor can therefore move to the new canonical maximizer of coordinate i. Every newly included factor vertex can also be visited while the other coordinates stay at their peaks. This proves the invariants and that every A_i is contained in the projection of C_t(x). Each strict expansion adds a factor vertex, so there are at most N-d strict expansions.

At a fixed point, suppose a product path from x first exits the product of the A_i. Its exiting step changes one coordinate i, from u_i in A_i to v_i outside A_i. Every other coordinate is still in its A_j and therefore has score at most p_j. Since the new product vertex is retained,

    f_i(v_i) >= t - sum_{j != i} p_j.

The old coordinate u_i belongs to the fixed factor component above that threshold. The edge u_i v_i therefore places v_i in the same component, a contradiction. No retained path leaves the product of the final sets. Together with the reachable-anchor argument, this proves that they are exactly the component projections.

Every reachable product vertex has coordinate scores bounded by p_i, while the joint anchor attains their sum. Thus it attains the component maximum.

Necessity of equal signatures follows from equality of connected components. For sufficiency, equal final coordinate sets have equal canonical coordinate maximizers. Both starting endpoints reach the same joint anchor. Reverse the path from y to that anchor, which is valid because every factor graph is undirected, and concatenate the paths. This proves connectivity. The result also makes the final sets independent of coordinate processing order. QED.

### Necessary warning about projections

Take two three-vertex paths, each scored (2,0,1), and t=2. From x=(0,0), both coordinate projections contain all three vertices. The retained vertex y=(2,2), scored 2, lies in that projection product but is isolated because each neighboring product vertex scores 1. Its own coordinate closure consists of two singletons. Thus neither a box representation of the component nor a one-sided projection membership test is valid.

## Proposition 2: a factor frontier algorithm

The same closure can be obtained by admitting one neighboring factor vertex at a time. Maintain a maximum-score heap of discovered but unadmitted neighbors for each A_i. Each factor vertex enters its heap at most once. Let q_i be the score of the top frontier vertex in factor i, let P=sum_i p_i, and let S=P-t. Store one representative per nonempty factor heap in a global minimum heap, with key

    k_i = p_i - q_i.

The top vertex in factor i is admissible precisely when k_i<=S. If any frontier vertex in any factor is admissible, a global minimum key is admissible. Admit that vertex, expose its unseen neighbors, and update only that factor's representative. If its score raises p_i by Delta, then S also rises by Delta. No other representative key changes.

Stop when every factor frontier is empty or the minimum key exceeds S. At that point no boundary vertex can enter any coordinate component at its final threshold, so each A_i is a fixed factor component. The invariant proof above applies to each vertex admission and gives the same exact projection closure.

After explicit input validation and ordering, this algorithm uses O((N+M) log(max(2,N))) arithmetic and comparison operations and O(N+M) storage, excluding an optional written path witness. Every admitted vertex exposes its adjacency list once, every factor vertex enters and leaves a local heap at most once, and every admission replaces at most one global representative. This bound is for a threshold query from one endpoint; connectivity requires two such queries. It does not enumerate all components or construct the full product merge tree.

The implementation retains its initially discovered parent for each frontier vertex. Because discovery occurs only from an admitted vertex, these parents form a spanning tree of each admitted coordinate set. Whenever its canonical peak changes, a tree path moves the old anchor coordinate to the new peak. All vertices on this path were admitted at thresholds at least as restrictive as the current threshold, so it is a valid product path when other coordinates are held at their peaks.

## Proposition 3: polynomial witness length

The anchor of factor i changes at most |V_i|-1 times, because each change accompanies admission of a new factor vertex. Each anchor change has a simple tree path of length at most |V_i|-1. Hence a path from x to its joint peak anchor uses at most

    sum_i (|V_i|-1)^2

product moves. Connected endpoints can therefore be joined by at most twice this quantity, by concatenating their paths through the common anchor. This is a constructive bound rather than a shortest-path claim. Paths can revisit vertices and can decrease the total score while remaining above t. The result does not assert monotone adaptation.

The present parent-tree reconstruction can walk toward the root even when the resulting path is short. It adds O(sum_i |V_i|^2) time, rather than a claimed bound proportional only to the number of emitted moves. Witness storage is bounded by O(sum_i (|V_i|-1)^2). Both bounds are polynomial in explicit factor input size even when the product has exponentially many vertices. Materializing each complete d-coordinate state adds a factor d to the size of that verbose output; the compact witness records coordinate moves.

## Corollary: exact rational widest capacity

Assume every input score is rational. Let D be the least common multiple of all score denominators and write w_i=D f_i, so every w_i is integral. First check whether x_i and y_i lie in the same baseline factor component for every i. If any do not, the product endpoints are disconnected before thresholding.

Otherwise let B_i be their common baseline factor component. The integer widest capacity lies in the inclusive interval

    L = sum_i min_{v in B_i} w_i(v),
    U = min(sum_i w_i(x_i), sum_i w_i(y_i)).

At L every vertex of the baseline product component is retained, so connectivity is certain. Connectivity is monotone in the threshold. Integer binary search therefore finds the greatest feasible threshold using O(log(1+U-L)) threshold queries. Divide the answer by D. The endpoint barrier is min(F(x),F(y)) minus this exact capacity. No product vertex enumeration or list of all attainable score sums is needed.

The numerical range must not be treated as a unit-cost input size. The bit length of D is at most the sum of the bit lengths of the input denominators. Scaled scores and their coordinate sums consequently have polynomial bit length in the rational input encoding. Additions, comparisons, least common multiples and binary search all have polynomial bit cost. The nearly linear operation count above concerns threshold decisions after normalization; it is not a claim of strongly polynomial bit complexity for arbitrary rational widest capacity computation.

## Established prior result and remaining novelty question

[Baryshnikov (2025), Appendix B, Proposition 7.1](https://link.springer.com/article/10.1007/s41468-025-00224-w) states that the merge tree of f(x)+g(y) depends only on the factor merge trees. The subsequent discussion describes the resulting bars and refers attachment rules to Le's earlier work. A [2022 author note](https://ymb.web.illinois.edu/2022/12/17/merge-tree-monoid/) states the same composition principle. Search-indexed material from [Le's 2023 Illinois thesis, Section 5.1](https://www.ideals.illinois.edu/items/127737/bitstreams/427840/data.pdf) explicitly describes construction for finite sums. The thesis landing page and PDF returned 403 during this audit, so its detailed attachment algorithm and asymptotic cost have not been checked.

Finite graphs fit the existing topological principle through piecewise linear edge extensions. On each product cell the sum is affine. Its retained corners connect to a maximizing corner along nondecreasing coordinate edges. Thus treating the input as a finite graph does not by itself establish a different structural contribution.

The possible additional statement is a direct factor-level endpoint query, its explicit operation and bit bounds, and a short constructive witness, without building the full composed merge tree. That algorithmic distinction is precise, but priority and novelty remain open until the closest construction is compared in full. The generic composition principle must be credited to prior work.

The distinction between one query and the full output matters. For d factors that are three-vertex paths scored (1,0,1), the threshold d has 2^d isolated retained vertices, each choosing one endpoint of every factor. Therefore an explicit list of components or all leaves of the composed superlevel merge tree can have exponential size in the number of factors. A polynomial query algorithm does not imply that the full object has a polynomial representation in every requested output format.

## Assumptions and scope

Undirected factors are needed when the second anchor path is reversed. Cartesian moves and additive scores are needed for the factor thresholds and the first-exit inequality. Coupled scores, simultaneous changes in several coordinates, directed transitions and unbounded or implicit factor graphs require separate arguments. Rational exactness in the implementation uses integers and Fraction inputs; floating point values are rejected. Disconnected factors, singleton factors, tied scores, equality at threshold and identical endpoints are included. The finite verification driver checks stated ranges only and does not substitute for this proof or for an independent proof audit.
