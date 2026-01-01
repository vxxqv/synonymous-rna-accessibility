"""Outcome-agnostic graph rules for the draft Cambray pilot.

No file readers or biological data access are provided by this module.
Unknown entries are None, never zero. Tolerances are explicit and default to0.
"""
import hashlib
import math
from fractions import Fraction


def finite(value):
    return value is not None and math.isfinite(value)


def _validate(graph, values):
    n = len(graph)
    if len(values) != n or not n or not values[0]:
        raise ValueError('Graph and trait dimensions disagree')
    traits = len(values[0])
    if any(len(row) != traits for row in values):
        raise ValueError('Ragged trait matrix')
    for v, neighbors in enumerate(graph):
        if len(set(neighbors)) != len(neighbors):
            raise ValueError('Duplicate neighbor')
        for w in neighbors:
            if not isinstance(w, int) or w < 0 or w >= n or w == v or v not in graph[w]:
                raise ValueError('Not a simple undirected graph')
    for row in values:
        for value in row:
            if value is not None and not finite(value):
                raise ValueError('Nonfinite known value')
    return traits


def reachable(graph, allowed, endpoints):
    s, t = endpoints
    if not (0 <= s < len(graph) and 0 <= t < len(graph)):
        raise ValueError('Endpoint out of range')
    allowed = set(allowed)
    if s not in allowed or t not in allowed:
        return False
    visited = {s}
    pending = [s]
    while pending:
        vertex = pending.pop()
        if vertex == t:
            return True
        for neighbor in graph[vertex]:
            if neighbor in allowed and neighbor not in visited:
                visited.add(neighbor)
                pending.append(neighbor)
    return False


def _thresholds(values, endpoints, thresholds):
    if len(endpoints) != 2 or any(not isinstance(v,int) or v < 0 or v >= len(values) for v in endpoints):
        raise ValueError('Invalid endpoints')
    s, t = endpoints
    if any(value is None for vertex in endpoints for value in values[vertex]):
        raise ValueError('Endpoint measurements must be known')
    if thresholds is None:
        return [min(a, b) for a, b in zip(values[s], values[t])]
    if len(thresholds) != len(values[0]) or any(not finite(v) for v in thresholds):
        raise ValueError('Invalid thresholds')
    return list(thresholds)


def _tolerances(tolerance, traits):
    tol = [tolerance] * traits if isinstance(tolerance, (int, float, Fraction)) else list(tolerance)
    if len(tol) != traits or any(not finite(v) or v < 0 for v in tol):
        raise ValueError('Invalid tolerance')
    return tol


def classify(graph, values, endpoints, thresholds=None, tolerance=0):
    traits = _validate(graph, values)
    if any(value is None for row in values for value in row):
        raise ValueError('Point classification requires complete values')
    tau = _thresholds(values, endpoints, thresholds)
    tol = _tolerances(tolerance, traits)
    allowed = [{v for v, row in enumerate(values) if row[k] >= tau[k] - tol[k]} for k in range(traits)]
    if any(v not in subset for v in endpoints for subset in allowed):
        return {'classification':'ineligible_endpoints','thresholds':tau,'individual_connectivity':None,'joint_connectivity':None}
    individual = [reachable(graph, subset, endpoints) for subset in allowed]
    joint = reachable(graph, set.intersection(*allowed), endpoints)
    label = 'jointly_connected' if joint else ('joint_only_obstruction' if all(individual) else 'single_trait_disconnection')
    return {'classification':label, 'thresholds':tau,'individual_connectivity':individual,'joint_connectivity':joint}


def interval_certificate(graph, intervals, endpoints, thresholds=None, tolerance=0):
    """Sufficient box-range certificates. None denotes an entirely unknown trait.

    The intervals are sets of possible values, not probability distributions.
    Unknown endpoints are refused. Complete but uncertain endpoints must be
    definitely acceptable; custom thresholds that violate this return ineligible.
    """
    lower, upper = [], []
    for row in intervals:
        lo, hi = [], []
        for interval in row:
            if interval is None:
                lo.append(None); hi.append(None)
            else:
                if len(interval) != 2 or not all(finite(v) for v in interval) or interval[0] > interval[1]:
                    raise ValueError('Invalid interval')
                lo.append(interval[0]); hi.append(interval[1])
        lower.append(lo); upper.append(hi)
    traits = _validate(graph, lower)
    _validate(graph, upper)
    tau = _thresholds(lower, endpoints, thresholds)
    tol = _tolerances(tolerance, traits)
    pessimistic = [{v for v, row in enumerate(lower) if row[k] is not None and row[k] >= tau[k]-tol[k]} for k in range(traits)]
    optimistic = [{v for v, row in enumerate(upper) if row[k] is None or row[k] >= tau[k]-tol[k]} for k in range(traits)]
    if any(v not in subset for v in endpoints for subset in pessimistic):
        return {'status':'ineligible_definite_endpoints','thresholds':tau}
    ip = [reachable(graph, subset, endpoints) for subset in pessimistic]
    io = [reachable(graph, subset, endpoints) for subset in optimistic]
    jp = reachable(graph, set.intersection(*pessimistic), endpoints)
    jo = reachable(graph, set.intersection(*optimistic), endpoints)
    if jp:
        status = 'certified_jointly_connected'
    elif not all(io):
        status = 'certified_single_trait_disconnection'
    elif all(ip) and not jo:
        status = 'certified_joint_only_obstruction'
    else:
        status = 'inconclusive'
    return {'status':status,'thresholds':tau,'individual_pessimistic':ip,'individual_optimistic':io,'joint_pessimistic':jp,'joint_optimistic':jo}


def missing_vertex_bounds(graph, values, endpoints, thresholds=None, tolerance=0):
    return interval_certificate(graph, [[None if v is None else (v,v) for v in row] for row in values], endpoints, thresholds, tolerance)


def component_tolerances(values, multiplier=1e-10):
    if not values or not values[0]:
        raise ValueError('Empty trait matrix')
    result = []
    for column in zip(*values):
        known = [abs(v) for v in column if v is not None]
        if not known or any(not finite(v) for v in known):
            raise ValueError('Trait has no finite values')
        result.append(multiplier * max(1, max(known)))
    return result


def additive_projection(values):
    """Corner order00,01,10,11; exact rational arithmetic is preserved."""
    if len(values) != 4 or any(not finite(v) for v in values):
        raise ValueError('Four finite square values required')
    contrast = values[0] - values[1] - values[2] + values[3]
    residual = [contrast/4, -contrast/4, -contrast/4, contrast/4]
    return {'contrast':contrast,'projection':[v-r for v,r in zip(values,residual)],'residual':residual}


def missing_corner_predictions(known, missing):
    """The withheld outcome must be absent, not passed then ignored."""
    if missing not in range(4) or set(known) != set(range(4)) - {missing}:
        raise ValueError('Require exactly the three known corners')
    if any(not finite(v) for v in known.values()):
        raise ValueError('Known corners must be finite')
    a, b, opposite = known[missing ^ 1], known[missing ^ 2], known[missing ^ 3]
    return {'additive':a+b-opposite,'adjacent_mean':(a+b)/2,'three_corner_mean':sum(known.values())/3}


def native_shortfalls(values, endpoints=(0,3)):
    if len(values) != 4 or any(not finite(v) for row in values for v in row):
        raise ValueError('Four complete finite corners required')
    tau = _thresholds(values, endpoints, None)
    return {v:[max(0,t-f) for t,f in zip(tau,row)] for v,row in enumerate(values) if v not in endpoints}


def threshold_geometry(graph, values, endpoints):
    """Exact critical atoms and open intervals over the endpoint-admissible domain.

    Values are converted to exact rational representations of their supplied
    numbers. No threshold-space probability or cell weight is assigned.
    """
    if _validate(graph, values) != 2 or any(v is None for row in values for v in row):
        raise ValueError('Complete two-trait field required')
    rational = [[Fraction(v) for v in row] for row in values]
    cap = _thresholds(rational,endpoints,None)
    axes = []
    for k in range(2):
        critical = sorted({row[k] for row in rational if row[k] <= cap[k]})
        cells = [{'kind':'open','lower':None,'upper':critical[0],'representative':critical[0]-1}]
        for i, value in enumerate(critical):
            cells.append({'kind':'atom','lower':value,'upper':value,'representative':value})
            if i+1 < len(critical):
                cells.append({'kind':'open','lower':value,'upper':critical[i+1],'representative':(value+critical[i+1])/2})
        axes.append(cells)
    cells = []
    for i, first in enumerate(axes[0]):
        for j, second in enumerate(axes[1]):
            result = classify(graph,rational,endpoints,[first['representative'],second['representative']])
            cells.append({'axis_indices':[i,j],'classification':result['classification']})
    return {'axes':axes,'cells':cells,'endpoint_caps':cap}


def corner_index(lexically_sorted_ids):
    if len(lexically_sorted_ids) != 4 or lexically_sorted_ids != sorted(set(lexically_sorted_ids)):
        raise ValueError('Four distinct lexically sorted IDs required')
    payload = 'C5-corner-v1|' + '|'.join(lexically_sorted_ids)
    return int.from_bytes(hashlib.sha256(payload.encode('utf-8')).digest(),'big') % 4
