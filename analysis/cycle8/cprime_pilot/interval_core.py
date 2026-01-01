"""Exact bounds for positive local-state proxies. No biological input access."""
from fractions import Fraction as F
import heapq

def interval(lo,hi):
    lo=F(lo);hi=None if hi is None else F(hi)
    if lo<0 or (hi is not None and (hi<=0 or hi<lo)):
        raise ValueError('Invalid positive-domain interval')
    return lo,hi

def widen(pair,gamma):
    lo,hi=interval(*pair);gamma=F(gamma)
    if gamma<1:raise ValueError('Widening cannot narrow')
    return lo/gamma,None if hi is None else hi*gamma

def marker_bounds(intervals):
    pairs=[interval(*p) for p in intervals]
    if not pairs:raise ValueError('Empty marker set')
    if len(pairs)==1:return {'minimum':F(1),'maximum':F(1),'maximum_attained':True}
    upper_finite=[u for l,u in pairs if u is not None]
    low=max(F(1),max(l for l,u in pairs)/min(upper_finite)) if upper_finite else F(1)
    high=F(1)
    for i,(_,u) in enumerate(pairs):
        for j,(l,_) in enumerate(pairs):
            if i==j:continue
            if u is None or l==0:return {'minimum':low,'maximum':None,'maximum_attained':False}
            high=max(high,u/l)
    return {'minimum':low,'maximum':high,'maximum_attained':True}

def common_order(contexts):
    """Exact feasibility of context-specific monotone maps of one free score."""
    if not contexts:raise ValueError('No contexts')
    vertices=sorted(next(iter(contexts.values())))
    if not vertices or any(sorted(d)!=vertices for d in contexts.values()):raise ValueError('Unequal state sets')
    fields={c:{v:interval(*d[v]) for v in vertices} for c,d in contexts.items()}
    edges={v:set() for v in vertices}; origins={}
    for c,d in sorted(fields.items()):
        for v in vertices:
            for w in vertices:
                if v!=w and d[v][1] is not None and d[v][1]<d[w][0]:
                    edges[v].add(w);origins.setdefault((v,w),[]).append({'context':c,'upper':d[v][1],'lower':d[w][0]})
    indegree={v:0 for v in vertices}
    for v in vertices:
        for w in edges[v]:indegree[w]+=1
    queue=[v for v in vertices if not indegree[v]];heapq.heapify(queue);order=[]
    while queue:
        v=heapq.heappop(queue);order.append(v)
        for w in sorted(edges[v]):
            indegree[w]-=1
            if not indegree[w]:heapq.heappush(queue,w)
    edge_records=[{'lower_state':v,'higher_state':w,'origins':origins[v,w]} for v,w in sorted(origins)]
    if len(order)!=len(vertices):
        color={v:0 for v in vertices};stack=[];cycle=None
        def dfs(v):
            nonlocal cycle
            color[v]=1;stack.append(v)
            for w in sorted(edges[v]):
                if color[w]==0:
                    if dfs(w):return True
                elif color[w]==1:
                    cycle=stack[stack.index(w):]+[w];return True
            stack.pop();color[v]=2;return False
        for v in vertices:
            if color[v]==0 and dfs(v):break
        if cycle is None:raise ArithmeticError('No directed-cycle witness')
        return {'status':'infeasible','edges':edge_records,'cycle':cycle,
                'cycle_edges':[{'lower_state':v,'higher_state':w,'origins':origins[v,w]} for v,w in zip(cycle,cycle[1:])]}
    positive=[x for d in fields.values() for p in d.values() for x in p if x is not None and x>0]
    delta=min(positive)/2 if positive else F(1)
    responses={}
    for c,d in fields.items():
        h=delta;responses[c]={}
        for v in order:
            l,u=d[v];h=max(h,l)
            if h<=0 or (u is not None and h>u):raise ArithmeticError('Invalid monotone completion')
            responses[c][v]=h
    return {'status':'feasible','edges':edge_records,'latent_order':order,'responses':responses}

def jsonable(obj):
    if isinstance(obj,F):return str(obj)
    if isinstance(obj,dict):return {k:jsonable(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):return [jsonable(x) for x in obj]
    return obj
