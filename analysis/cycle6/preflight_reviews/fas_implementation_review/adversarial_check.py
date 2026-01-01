"""Synthetic-only independent ordinal and threshold-connectivity oracles.

No biological data file is accepted. Candidate module lives in a review copy.
"""
import argparse
import hashlib
import importlib.util
import itertools as it
import json
import math
from fractions import Fraction as F
from pathlib import Path
from collections import deque, Counter


def require(value, message):
    if not value:
        raise RuntimeError(message)


def cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def convex_origin(rows):
    """Exact origin-in-convex-hull oracle in R2 by Caratheodory supports <=3."""
    rows = sorted(set(rows))
    if (0,0) in rows:
        return True
    for a,b in it.combinations(rows,2):
        if cross(a,b)==0 and sum(u*v for u,v in zip(a,b)) < 0:
            return True
    for a,b,c in it.combinations(rows,3):
        det=cross(tuple(v-u for u,v in zip(a,b)),tuple(v-u for u,v in zip(a,c)))
        if det:
            weights=[F(cross(b,c),det),F(cross(c,a),det),F(cross(a,b),det)]
            require(sum(weights)==1,'Barycentric identity')
            if min(weights)>=0:
                return True
    return False


def reference_status(x, lower, upper):
    # All finite strict inequalities can be scaled simultaneously to unit margin.
    rows=[(x[j][-2]-x[i][-2],x[j][-1]-x[i][-1])
          for i in range(len(x)) for j in range(len(x))
          if upper[i] is not None and upper[i]<lower[j]]
    return 'infeasible' if convex_origin(rows) else 'feasible'


def check_certificate(x,lower,upper,r):
    n=len(x)
    pairs={(i,j) for i in range(n) for j in range(n)
           if upper[i] is not None and upper[i]<lower[j]}
    require(r['forced_orders']==len(pairs),'Wrong order count')
    require(r['exact_checked'] is True,'Unchecked return')
    if r['status']=='feasible':
        b=list(map(F,r['beta'])); latent=list(map(F,r['latent'])); y=list(map(F,r['response']))
        require(len(b)==len(x[0]) and len(latent)==n and len(y)==n,'Wrong certificate dimension')
        require(latent==[sum(a*z for a,z in zip(row,b)) for row in x],'Invalid latent witness')
        require(all(lower[i]<=y[i] and (upper[i] is None or y[i]<=upper[i]) for i in range(n)),'Outside intervals')
        require(all(y[i]<=y[j] for i in range(n) for j in range(n) if latent[i]<=latent[j]),'Not a single-valued monotone map')
        require(all(latent[i]<latent[j] for i,j in pairs),'Bad strict order')
    elif r['status']=='infeasible':
        total=F(0); cancel=[F(0)]*len(x[0])
        for d in r['dual']:
            i,j=d['lower_vertex'],d['higher_vertex'];w=F(d['weight'])
            require((i,j) in pairs and w>=0,'Dual includes invalid order/weight')
            delta=[v-u for u,v in zip(x[i],x[j])]
            require(d['difference']==delta,'Dual difference not independently verifiable')
            total+=w;cancel=[a+w*z for a,z in zip(cancel,delta)]
        require(total==1 and all(z==0 for z in cancel),'Dual is not Farkas cancellation')
    else:
        raise RuntimeError('Unknown status')


def reachable(scores,edges,start,target,t):
    if scores[start]<t or scores[target]<t:return False
    seen={start};queue=deque([start])
    while queue:
        v=queue.popleft()
        for a,b in edges:
            nxt=b if a==v else a if b==v else None
            if nxt is not None and nxt not in seen and scores[nxt]>=t:
                seen.add(nxt);queue.append(nxt)
    return target in seen


def capacity(scores,edges,a,b):
    available=[t for t in set(scores) if reachable(scores,edges,a,b,t)]
    return max(available) if available else -math.inf


def extended(s):
    return float(s) if s in ('inf','-inf') else F(s)


def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();require(not args.out.exists(),'No overwrite')
    raw=args.candidate.read_bytes();sha=hashlib.sha256(raw).hexdigest()
    spec=importlib.util.spec_from_file_location('held_target',args.candidate)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    counts=Counter();failures=[]
    matrices=[[(1,0,0),(1,1,0),(1,2,0)],[(1,0,0),(1,1,0),(1,0,1)],
              [(1,0,0),(1,0,0),(1,1,0)],[(1,0,0),(1,1,1),(1,-1,-1)]]
    grid=[(F(0),F(0)),(F(0),F(1)),(F(1),F(1)),(F(1),F(2)),(F(0),None),(F(2),None)]
    for x in matrices:
        for intervals in it.product(grid,repeat=3):
            lower=[a for a,b in intervals];upper=[b for a,b in intervals]
            r=mod.order_certificate(x,lower,upper)
            require(r['status']==reference_status(x,lower,upper),'Convex-hull oracle mismatch')
            check_certificate(x,lower,upper,r);counts['ordinal_cases']+=1;counts[r['status']]+=1
    q2=list(it.product(range(2),repeat=2));x2=[[1,a,b] for a,b in q2]
    for values in it.product((F(-2),F(0),F(3)),repeat=4):
        for perm in [tuple(range(4)),(3,0,2,1),(2,3,1,0)]:
            x=[x2[i] for i in perm];v=[values[i] for i in perm]
            r=mod.order_certificate(x,v,v)
            require(r['status']==reference_status(x,v,v),'Q2 order mismatch')
            check_certificate(x,v,v,r);counts['ordinal_cases']+=1;counts[r['status']]+=1
    # Full declared factor shape, with constructed exact additive latent and arbitrary nondecreasing plateau.
    states=list(it.product(range(3),range(3),range(2),range(2),range(2)))
    x=mod.features(states)
    beta=[0,2,-3,1,4,2,-1,3]
    latent=[sum(a*b for a,b in zip(row,beta)) for row in x]
    vals=[F(max(0,t)**3,7) for t in latent]
    for perm in [list(range(72)),list(reversed(range(72))),list(range(0,72,2))+list(range(1,72,2))]:
        xx=[x[i] for i in perm];vv=[vals[i] for i in perm]
        r=mod.order_certificate(xx,vv,vv);require(r['status']=='feasible','Constructed monotone additive field rejected')
        check_certificate(xx,vv,vv,r);counts['full72_feasible']+=1
    # 64 labelled graphs, including empty/disconnected/self-loop/duplicate-edge variants.
    potential=list(it.combinations(range(4),2))
    graph_list=[[e for k,e in enumerate(potential) if mask & (1<<k)] for mask in range(64)]
    for edges in graph_list:
        for vals in it.product((F(-1),F(1),F(2)),repeat=4):
            d=mod.widest_paths(vals,edges)
            for a,b in it.combinations(range(4),2):
                require(d[a][b]==capacity(vals,edges,a,b),'False widest path')
                counts['widest_endpoint_checks']+=1
    interval_grid=[(F(0),F(0)),(F(0),F(2)),(F(1),F(2)),(F(1),None)]
    special_graphs=[[],[(0,1),(1,3),(0,2),(2,3)],[(0,1),(1,2),(2,3)],[(0,1),(2,3)]]
    for intervals in it.product(interval_grid,repeat=4):
        lower=[a for a,b in intervals];upper=[b for a,b in intervals];upper_scores=[math.inf if b is None else b for b in upper]
        for edges in special_graphs:
            r=mod.bounded_pairs(lower,upper,edges)
            for row in r:
                a,b=row['a'],row['b'];cl=capacity(lower,edges,a,b);cu=capacity(upper_scores,edges,a,b)
                endpoint=min(lower[a],lower[b]);width=max(F(0),endpoint-cu)
                require(extended(row['bottleneck_lower'])==cl and extended(row['bottleneck_upper'])==cu,'Wrong interval capacity')
                require(extended(row['certified_interval_width'])==width,'Wrong certified width')
                require(row['certified_nonempty']==(cu<endpoint),'Wrong open/closed interval status')
                if cu==endpoint:require(not row['certified_nonempty'],'Tie falsely blocks')
                # A finite midpoint in a positive certified threshold interval must disconnect the optimistic field.
                if math.isfinite(cu) and cu<endpoint:
                    t=(cu+endpoint)/2
                    require(not reachable(upper_scores,edges,a,b,t),'False optimistic disconnection')
                counts['interval_endpoint_checks']+=1
    # Boundary risks are retained as findings, never suppressed as failed intended examples.
    edges=[(0,1),(1,2)];l=[2,0,2];u=[2,3,2]
    expected=mod.bounded_pairs(l,u,edges)
    try:
        actual=mod.bounded_pairs(l,u,iter(edges))
        if actual!=expected:failures.append({'case':'edge_iterator_exhaustion','expected_02':expected[1],'actual_02':actual[1]})
    except (TypeError,ValueError):
        counts['explicit_iterator_refusal']+=1
    invalid_cases=[('noninteger_features',lambda:mod.order_certificate([[0.1],[0.2]],[0,1],[0,1])),
                   ('invalid_interval',lambda:mod.bounded_pairs([2,0],[1,1],[(0,1)])),
                   ('invalid_edge',lambda:mod.bounded_pairs([0,0],[1,1],[(0,2)])),
                   ('nan_score',lambda:mod.widest_paths([0,float('nan')],[(0,1)]))]
    for name,call in invalid_cases:
        try:call()
        except (ValueError,TypeError,OverflowError):counts['invalid_inputs_refused']+=1
        else:failures.append({'case':name,'error':'Invalid input accepted'})
    require(raw==args.candidate.read_bytes(),'Candidate source changed')
    payload={'status':'BLOCKED' if failures else 'PASS','candidate_sha256':sha,
             'counts':dict(counts),'failures':failures,
             'scope':'Synthetic only. Exact convex-hull ordinal oracle and threshold BFS; no author numeric reference imported.',
             'certificate_limits':'Fail-closed ArithmeticError remains allowed for unresolvable floating discovery; every returned certificate checked rationally.',
             'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.out.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(payload,indent=2))


if __name__=='__main__':
    main()
