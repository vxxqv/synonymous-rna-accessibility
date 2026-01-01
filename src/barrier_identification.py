"""Endpoint barrier identification on finite undirected graphs.

Intervals are deterministic simultaneous bounds, not marginal confidence
intervals. Sharpness requires every vector in their Cartesian product to be
admissible. Widest-path search is established; the endpoint-separated box
calculation is the object implemented here.
"""
from heapq import heappush, heappop
from math import inf, isfinite


def _validate(adjacency, values, a, b):
    n=len(adjacency)
    if len(values)!=n or not 0<=a<n or not 0<=b<n:
        raise ValueError('Graph, score and endpoint dimensions disagree')
    if not all(isfinite(x) for x in values):
        raise ValueError('Scores and interval bounds must be finite')
    for i, neighbors in enumerate(adjacency):
        if any(not 0<=j<n or i not in adjacency[j] for j in neighbors):
            raise ValueError('Expected a finite undirected adjacency list')


def _widest(adjacency, values, a, b):
    capacity=[-inf]*len(adjacency);capacity[a]=values[a]
    heap=[(-capacity[a],a)]
    while heap:
        neg,u=heappop(heap)
        if -neg!=capacity[u]:continue
        if u==b:return capacity[u]
        for v in adjacency[u]:
            new=min(capacity[u],values[v])
            if new>capacity[v]:
                capacity[v]=new;heappush(heap,(-new,v))
    raise ValueError('Endpoints are in different baseline components')


def bottleneck(adjacency, values, a, b):
    _validate(adjacency,values,a,b)
    return _widest(adjacency,values,a,b)


def barrier(adjacency, values, a, b):
    capacity=bottleneck(adjacency,values,a,b)
    return min(values[a],values[b])-capacity


def internal_bottleneck(adjacency, values, a, b):
    _validate(adjacency,values,a,b)
    if a==b:return inf
    modified=list(values);modified[a]=modified[b]=inf
    return _widest(adjacency,modified,a,b)


def sharp_barrier_bounds(adjacency, lower, upper, a, b):
    _validate(adjacency,lower,a,b);_validate(adjacency,upper,a,b)
    if any(l>u for l,u in zip(lower,upper)):
        raise ValueError('Lower bound exceeds upper bound')
    if a==b:return (0,0)
    low=max(0,min(lower[a],lower[b])-internal_bottleneck(adjacency,upper,a,b))
    high=max(0,min(upper[a],upper[b])-internal_bottleneck(adjacency,lower,a,b))
    return low,high


def threshold_status(adjacency, lower, upper, a, b, threshold):
    """Classify a fixed threshold, or abstain; endpoints must be retained."""
    _validate(adjacency,lower,a,b);_validate(adjacency,upper,a,b)
    if not isfinite(threshold):
        raise ValueError('Threshold must be finite')
    if any(l>u for l,u in zip(lower,upper)):
        raise ValueError('Lower bound exceeds upper bound')
    # Baseline connectivity is required even when one endpoint is uncertain.
    _widest(adjacency,[0]*len(lower),a,b)
    if min(lower[a],lower[b])<threshold:return 'endpoints_not_certified'
    if _widest(adjacency,lower,a,b)>=threshold:return 'connected'
    if _widest(adjacency,upper,a,b)<threshold:return 'separated'
    return 'unresolved'


def additive_background_barrier(local_barrier, initial_background, maximum_background):
    """Requires a superlevel path from the initial background to its maximum."""
    if not all(isfinite(x) for x in (local_barrier,initial_background,maximum_background)):
        raise ValueError('Barrier and background scores must be finite')
    if local_barrier<0 or maximum_background<initial_background:
        raise ValueError('Invalid barrier or background ordering')
    return max(0,local_barrier-(maximum_background-initial_background))
