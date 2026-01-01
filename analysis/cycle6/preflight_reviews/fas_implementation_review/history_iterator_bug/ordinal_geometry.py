"""Exact-checked ordinal null and bounded accessibility; no data acquisition."""
from fractions import Fraction as F
import itertools
import math
import numpy as np
from scipy.optimize import linprog


def features(states):
    """Treatment-coded categorical main effects, including an intercept."""
    levels = [sorted(set(s[k] for s in states)) for k in range(len(states[0]))]
    return [[1] + [int(s[k] == x) for k, a in enumerate(levels) for x in a[1:]]
            for s in states]


def order_certificate(x, lower, upper):
    """None upper means positive infinity. All asserted certificates are rational."""
    if not x or not x[0] or any(z != int(z) for row in x for z in row):
        raise ValueError('A nonempty integer design matrix is required')
    x = [[int(z) for z in row] for row in x]
    lower = [F(z) for z in lower]
    upper = [None if z is None else F(z) for z in upper]
    n = len(x)
    if len(lower) != n or len(upper) != n or n == 0:
        raise ValueError('Empty or inconsistent observations')
    if len({len(r) for r in x}) != 1:
        raise ValueError('Inconsistent feature dimensions')
    if any(u is not None and l > u for l, u in zip(lower, upper)):
        raise ValueError('Empty interval')
    pairs = [(i,j) for i in range(n) for j in range(n)
             if upper[i] is not None and upper[i] < lower[j]]
    a = [[v-u for u,v in zip(x[i],x[j])] for i,j in pairs]
    if not a:
        beta = [F(0)] * len(x[0])
        result = None
    else:
        result = linprog(np.zeros(len(x[0])), A_ub=-np.array(a,dtype=float),
                         b_ub=-np.ones(len(a)), bounds=[(None,None)]*len(x[0]), method='highs')
        if result.status == 0:
            beta = [F(float(v)).limit_denominator(1000000) for v in result.x]
            margins = [sum(v*b for v,b in zip(row,beta)) for row in a]
            if min(margins) <= 0:
                raise ArithmeticError('Floating feasible result lacks exact strict certificate')
            scale = min(margins)
            beta = [v/scale for v in beta]
        elif result.status != 2:
            raise ArithmeticError(f'Unresolved optimizer status {result.status}')
    if result is None or result.status == 0:
        latent = [sum(v*b for v,b in zip(row,beta)) for row in x]
        y = [max(lower[j] for j in range(n) if latent[j] <= latent[i]) for i in range(n)]
        if not all(upper[i] is None or y[i] <= upper[i] for i in range(n)):
            raise ArithmeticError('Exact monotone-response certificate failed')
        if not all(latent[i] < latent[j] for i,j in pairs):
            raise ArithmeticError('Exact forced-order certificate failed')
        return {'status':'feasible','forced_orders':len(pairs),'beta':list(map(str,beta)),
                'latent':list(map(str,latent)),'response':list(map(str,y)),
                'exact_checked':True}
    # A rational Farkas witness is a convex combination of difference rows equal to zero.
    ae = np.vstack([np.array(a,dtype=float).T, np.ones(len(a))])
    rhs = np.r_[np.zeros(len(x[0])),1.]
    dual = linprog(np.zeros(len(a)), A_eq=ae, b_eq=rhs, bounds=(0,None), method='highs')
    if not dual.success:
        raise ArithmeticError('No dual certificate found')
    support = [k for k,v in enumerate(dual.x) if v > 1e-9]
    # Solve only the sparse support exactly, independently of floating residuals.
    from sympy import Matrix
    mat = Matrix([[a[k][j] for k in support] for j in range(len(x[0]))]
                 + [[1]*len(support)])
    sol, params = mat.gauss_jordan_solve(Matrix([0]*len(x[0])+[1]))
    if params.rows:
        sol = sol.subs({v:0 for v in params})
    weights = [F(int(v.p),int(v.q)) for v in sol]
    if not all(v >= 0 for v in weights) or sum(weights) != 1:
        raise ArithmeticError('Invalid exact dual weights')
    if any(sum(w*a[k][j] for k,w in zip(support,weights)) != 0 for j in range(len(x[0]))):
        raise ArithmeticError('Exact cancellation failed')
    return {'status':'infeasible','forced_orders':len(pairs),'exact_checked':True,
            'dual':[{'lower_vertex':pairs[k][0],'higher_vertex':pairs[k][1],
                     'weight':str(w),'difference':a[k]} for k,w in zip(support,weights) if w]}


def widest_paths(scores, edges):
    """Floyd max-min closure. Unknown upper scores use math.inf."""
    n=len(scores)
    if n == 0 or any(math.isnan(float(s)) for s in scores):
        raise ValueError('Nonempty scores without NaN required')
    if any(i != int(i) or j != int(j) or not (0 <= i < n and 0 <= j < n) for i,j in edges):
        raise ValueError('Invalid edge index')
    d=[[-math.inf]*n for _ in range(n)]
    for i,s in enumerate(scores): d[i][i]=s
    for i,j in edges: d[i][j]=d[j][i]=min(scores[i],scores[j])
    for k in range(n):
        for i in range(n):
            for j in range(n):
                d[i][j]=max(d[i][j], min(d[i][k],d[k][j]))
    return d


def bounded_pairs(lower, upper, edges):
    if len(lower) != len(upper) or not lower:
        raise ValueError('Nonempty equal-length intervals required')
    lower=[F(z) for z in lower]
    upper=[None if z is None else F(z) for z in upper]
    if any(u is not None and l > u for l,u in zip(lower,upper)):
        raise ValueError('Empty interval')
    lo=widest_paths(lower,edges)
    hi=widest_paths([math.inf if z is None else z for z in upper],edges)
    return [{'a':i,'b':j,'endpoint_lower':str(min(lower[i],lower[j])),
             'bottleneck_lower':str(lo[i][j]),'bottleneck_upper':str(hi[i][j]),
             'certified_interval_width':str(max(F(0),min(lower[i],lower[j])-hi[i][j])),
             'certified_nonempty':hi[i][j] < min(lower[i],lower[j]),
             'threshold_convention':'(bottleneck_upper, endpoint_lower] if nonempty; acceptance >= threshold'}
            for i,j in itertools.combinations(range(len(lower)),2)]


def graph_edges(states, factor_edges):
    ans=[]
    for i,j in itertools.combinations(range(len(states)),2):
        diff=[k for k,(a,b) in enumerate(zip(states[i],states[j])) if a!=b]
        if len(diff)==1:
            k=diff[0]
            if frozenset((states[i][k],states[j][k])) in factor_edges[k]: ans.append((i,j))
    return ans
