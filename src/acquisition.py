"""Exact vertex queries for threshold connectivity on finite undirected graphs.

Query costs count newly observed internal vertices. Endpoints are known present.
The policies are references, not claimed optimal. Unknown is a separate state
and is never interpreted as absence.
"""
from collections import deque
from functools import lru_cache
from heapq import heappop, heappush

import networkx as nx


POLICIES = ('random', 'degree', 'lazy_forward', 'lazy_bisection', 'path_cut')


def validate(adjacency, state, a, b):
    n = len(adjacency)
    if len(state) != n or not 0 <= a < n or not 0 <= b < n:
        raise ValueError('Graph, state and endpoint dimensions disagree')
    if any(value not in (-1, 0, 1) for value in state):
        raise ValueError('States must be unknown (-1), absent (0), or present (1)')
    if state[a] != 1 or state[b] != 1:
        raise ValueError('Both endpoints must be known present')
    for u, neighbors in enumerate(adjacency):
        if len(set(neighbors)) != len(neighbors) or u in neighbors:
            raise ValueError('Expected a simple graph')
        if any(not 0 <= v < n or u not in adjacency[v] for v in neighbors):
            raise ValueError('Expected a finite undirected graph')


def shortest_path(adjacency, allowed, a, b):
    """Lexicographic BFS path, or None when no allowed path exists."""
    if not allowed[a] or not allowed[b]:
        return None
    parents = {a: None}
    queue = deque([a])
    while queue:
        u = queue.popleft()
        if u == b:
            path = []
            while u is not None:
                path.append(u)
                u = parents[u]
            return path[::-1]
        for v in sorted(adjacency[u]):
            if allowed[v] and v not in parents:
                parents[v] = u
                queue.append(v)
    return None


def certificate(adjacency, state, a, b):
    """Return an exact conditional certificate, or unresolved."""
    validate(adjacency, state, a, b)
    if shortest_path(adjacency, [s == 1 for s in state], a, b) is not None:
        return 'connected'
    if shortest_path(adjacency, [s != 0 for s in state], a, b) is None:
        return 'separated'
    return 'unresolved'


def _vertex_cut(adjacency, allowed, costs, a, b, prohibited):
    """Node splitting with integral capacities and sorted insertion order."""
    graph = nx.DiGraph()
    for u in range(len(adjacency)):
        if allowed[u]:
            graph.add_edge(2*u, 2*u+1, capacity=costs[u])
    for u in range(len(adjacency)):
        if allowed[u]:
            for v in sorted(adjacency[u]):
                if allowed[v]:
                    graph.add_edge(2*u+1, 2*v, capacity=prohibited)
    value, (reachable, unreachable) = nx.minimum_cut(
        graph, 2*a+1, 2*b, capacity='capacity',
        flow_func=nx.algorithms.flow.edmonds_karp)
    cut = [u for u in range(len(adjacency))
           if 2*u in reachable and 2*u+1 in unreachable]
    return value, cut


def minimum_unknown_separator(adjacency, state, a, b):
    """Minimum cardinality separator, with lexicographically first tied set.

    Integral weight 2**n - 2**(n-1-v) gives cardinality strict priority and
    uniquely favors the lower indexed member of the first differing pair.
    Known present vertices and graph edges have prohibitive capacity.
    """
    validate(adjacency, state, a, b)
    if certificate(adjacency, state, a, b) != 'unresolved':
        raise ValueError('A query is needed only in an unresolved state')
    n = len(state)
    unit = 1 << n
    prohibited = n * unit + 1
    costs = [unit - (1 << (n-1-v)) if s == -1 else prohibited
             for v, s in enumerate(state)]
    value, cut = _vertex_cut(adjacency, [s != 0 for s in state],
                             costs, a, b, prohibited)
    if value >= prohibited or not cut or any(state[v] != -1 for v in cut):
        raise AssertionError('No valid unknown separator was recovered')
    removed = set(cut)
    if shortest_path(adjacency, [s != 0 and v not in removed
                                for v, s in enumerate(state)], a, b):
        raise AssertionError('Recovered separator does not disconnect endpoints')
    return cut


def choose_query(adjacency, state, a, b, policy, step, rng):
    if policy not in POLICIES:
        raise ValueError('Unknown policy')
    if certificate(adjacency, state, a, b) != 'unresolved':
        raise ValueError('An additional query is unnecessary')
    unknown = [v for v, s in enumerate(state) if s == -1]
    if policy == 'random':
        return int(rng.choice(unknown))
    if policy == 'degree':
        return min(unknown, key=lambda v: (-len(adjacency[v]), v))
    if policy == 'path_cut' and step % 2 == 1:
        return min(minimum_unknown_separator(adjacency, state, a, b))
    path = shortest_path(adjacency, [s != 0 for s in state], a, b)
    candidates = [v for v in path if state[v] == -1]
    if policy == 'lazy_forward':
        return candidates[0]
    # Central unqueried position on the complete possible path. Index, not
    # path order, resolves equal distance from the path midpoint.
    locations = {v: i for i, v in enumerate(path)}
    return min(candidates, key=lambda v: (abs(2*locations[v]-(len(path)-1)), v))


def clairvoyant_additional_cost(adjacency, state, truth, a, b):
    """Minimum true certificate size after initial observations.

    This lower bound uses the complete truth and is not a query policy.
    """
    validate(adjacency, state, a, b)
    if any(t not in (0, 1) for t in truth) or len(truth) != len(state):
        raise ValueError('Truth must be a complete binary state')
    if any(s != -1 and s != t for s, t in zip(state, truth)):
        raise ValueError('Observed state disagrees with truth')
    if certificate(adjacency, state, a, b) != 'unresolved':
        return 0
    if shortest_path(adjacency, truth, a, b) is not None:
        costs = [int(s == -1) for s in state]
        distance = {a: 0}
        queue = [(0, a)]
        while queue:
            cost, u = heappop(queue)
            if cost != distance[u]:
                continue
            if u == b:
                return cost
            for v in adjacency[u]:
                if truth[v]:
                    candidate = cost + costs[v]
                    if candidate < distance.get(v, len(state)+1):
                        distance[v] = candidate
                        heappush(queue, (candidate, v))
        raise AssertionError('Connected truth has no path')
    prohibited = len(state) + 1
    costs = [1 if t == 0 and s == -1 else prohibited
             for t, s in zip(truth, state)]
    value, _ = _vertex_cut(adjacency, [s != 0 for s in state],
                            costs, a, b, prohibited)
    if value >= prohibited:
        raise AssertionError('Separated truth has no finite negative certificate')
    return int(value)


def exact_minimax(adjacency, a, b):
    """Return a cached exact worst case recurrence for small query graphs.

    C(s)=0 at certificates, otherwise 1+min_v max(C(s[v=0]),C(s[v=1])).
    Callers must bound graph size before using the exponential recurrence.
    """
    if len(adjacency) > 10:
        raise ValueError('Exact minimax reference is capped at ten vertices')

    @lru_cache(None)
    def solve(state):
        if certificate(adjacency, state, a, b) != 'unresolved':
            return 0
        return 1 + min(max(solve(state[:v]+(0,)+state[v+1:]),
                           solve(state[:v]+(1,)+state[v+1:]))
                       for v, s in enumerate(state) if s == -1)
    return solve
