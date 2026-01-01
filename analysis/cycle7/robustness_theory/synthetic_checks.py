"""Synthetic exact checks only; never loads an observed landscape.

Rational Fourier elimination for tiny fixed examples, independently compared
with the three possible latent orders of a one-dimensional rank-one design.
"""
import sys
if not __debug__:
    raise RuntimeError('Optimized execution refused')
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json


def check(ok, message):
    if not ok:
        raise AssertionError(message)


def linear(rows, p):
    if p == 0:
        return () if all(rhs <= 0 for _,rhs in rows) else None
    positive=[(a,r) for a,r in rows if a[-1]>0]
    negative=[(a,r) for a,r in rows if a[-1]<0]
    smaller=[(a[:-1],r) for a,r in rows if a[-1]==0]
    for a,r in positive:
        for b,s in negative:
            smaller.append((tuple(a[k]/a[-1]-b[k]/b[-1] for k in range(p-1)),r/a[-1]-s/b[-1]))
    x=linear(smaller,p-1)
    if x is None:
        return None
    lows=[(r-sum(z*v for z,v in zip(a[:-1],x)))/a[-1] for a,r in positive]
    highs=[(r-sum(z*v for z,v in zip(a[:-1],x)))/a[-1] for a,r in negative]
    v=max(lows) if lows else (min(highs) if highs else Q(0))
    if highs and v>min(highs):
        return None
    out=x+(v,)
    check(all(sum(z*v for z,v in zip(a,out))>=r for a,r in rows),'Recovered rational witness failed')
    return out


def valid(L,U,a):
    check(len(L)==len(U)==len(a)>0,'Nonempty matching dimensions')
    check(all(type(z) is int and z in (0,1) for z in a),'Weights must be binary')
    check(all(l>=0 and (u is None or u>=l) for l,u in zip(L,U)),'Nonnegative valid intervals required')


def bounds(L,U,a,metric,z):
    """Monomial (c,p) denotes c*sqrt(z)^p for mult; rational (c,0) for add."""
    if metric=='add':
        return [(max(Q(0),l-w*z),0) for l,w in zip(L,a)], [None if u is None else (u+w*z,0) for u,w in zip(U,a)]
    return [((l,0) if w==0 else (l/z,1)) for l,w in zip(L,a)], [None if u is None else ((u,0) if w==0 else (u,1)) for u,w in zip(U,a)]


def key(x,metric,z):
    c,p=x
    check(c>=0,'Negative monomial cannot be compared by squaring')
    return c*c*z**p if metric=='mult' else c


def active(L,U,a,metric,z):
    lower,upper=bounds(L,U,a,metric,z)
    return [(i,j) for i in range(len(L)) for j in range(len(L))
            if upper[i] is not None and key(upper[i],metric,z)<key(lower[j],metric,z)]


def events(L,U,a,metric):
    valid(L,U,a)
    candidates={Q(0) if metric=='add' else Q(1)}
    for i,j in product(range(len(L)),repeat=2):
        k=a[i]+a[j]
        if U[i] is None or k==0 or L[j]<=U[i]:
            continue
        if metric=='add':
            candidates.add((L[j]-U[i])/k)
        elif U[i]>0:
            candidates.add((L[j]/U[i])**(2//k))
    return sorted(candidates)


def fit(X,L,U,a,metric,z):
    pairs=active(L,U,a,metric,z)
    rows=[(tuple(Q(b-c) for b,c in zip(X[j],X[i])),Q(1)) for i,j in pairs]
    beta=linear(rows,len(X[0]))
    if beta is None:
        return None
    latent=[sum(Q(x)*v for x,v in zip(row,beta)) for row in X]
    lower,upper=bounds(L,U,a,metric,z)
    y=[max((lower[j] for j in range(len(L)) if latent[j]<=latent[i]),key=lambda v:key(v,metric,z)) for i in range(len(L))]
    check(all(key(lower[i],metric,z)<=key(y[i],metric,z) and (upper[i] is None or key(y[i],metric,z)<=key(upper[i],metric,z)) for i in range(len(L))),'Exact interval lift')
    check(all(key(y[i],metric,z)<=key(y[j],metric,z) for i in range(len(L)) for j in range(len(L)) if latent[i]<=latent[j]),'Monotonicity including ties')
    return {'beta':list(map(str,beta)),'latent':list(map(str,latent)),
            'response':[{'coefficient':str(c),'sqrt_t_power':p} for c,p in y]}


def threshold(X,L,U,a,metric):
    for z in events(L,U,a,metric):
        answer=fit(X,L,U,a,metric,z)
        if answer is not None:
            return z,answer
    return None,None


def order_oracle(L,U,a,metric,z):
    # X=(0),(1),(2): only increasing, decreasing, and all-tied latent orders.
    lower,upper=bounds(L,U,a,metric,z)
    for latent in [(0,1,2),(2,1,0),(0,0,0)]:
        if all(upper[i] is None or max(key(lower[j],metric,z) for j in range(3) if latent[j]<=latent[i])<=key(upper[i],metric,z) for i in range(3)):
            return True
    return False


def main():
    X=[(Q(0),),(Q(1),),(Q(2),)]
    cases=[]
    fixtures=[
      ('peak_unanchored',[1,4,2],[1,4,2],[1,1,1],Q(1),Q(2)),
      ('peak_middle_anchored',[1,4,2],[1,4,2],[1,0,1],Q(2),Q(4)),
      ('constant_compatible',[1,1,1],[1,1,1],[1,0,1],Q(0),Q(1)),
      ('zero_persistent',[0,1,0],[0,1,0],[1,1,1],Q(1,2),None),
      ('unknown_all',[0,0,0],[None,None,None],[1,0,1],Q(0),Q(1)),
    ]
    for name,l,u,a,ea,em in fixtures:
        L=list(map(Q,l));U=[None if v is None else Q(v) for v in u]
        za,ba=threshold(X,L,U,a,'add');zm,bm=threshold(X,L,U,a,'mult')
        check(za==ea and zm==em,'Known answer '+name)
        cases.append({'name':name,'additive_epsilon':None if za is None else str(za),'multiplicative_t':None if zm is None else str(zm),'additive_boundary_witness':ba,'multiplicative_boundary_witness':bm})
    # Duplicate-feature contradictions: fragile first certificate is not global radius.
    Y=[(Q(0),),(Q(0),),(Q(1),),(Q(1),)]
    L=list(map(Q,[1,2,1,4]));a=[1]*4
    za,_=threshold(Y,L,L,a,'add');zm,_=threshold(Y,L,L,a,'mult')
    check(za==Q(3,2) and zm==4,'Full-system versus first witness')
    # Permanent pair between fixed duplicate rows precludes every finite widening.
    for metric in ('add','mult'):
        z,_=threshold([(Q(0),)]*2,[Q(1),Q(2)],[Q(1),Q(2)],[0,0],metric)
        check(z is None,'Fixed incompatible anchors must stay infeasible')
    comparisons=0
    squared_comparisons=0
    for vals in product([Q(0),Q(1),Q(2),Q(4)],repeat=3):
      for a in product([0,1],repeat=3):
       for metric in ('add','mult'):
        for z in events(vals,vals,a,metric):
            check((fit(X,vals,vals,a,metric,z) is not None)==order_oracle(vals,vals,a,metric,z),'Independent weak-order disagreement')
            comparisons+=1
            if metric=='mult':
                squared_L=[v*v/z**w for v,w in zip(vals,a)]
                squared_U=[v*v*z**w for v,w in zip(vals,a)]
                check(active(vals,vals,a,'mult',z)==active(squared_L,squared_U,[0]*3,'add',Q(0)),'Squared order correspondence')
                first=fit(X,vals,vals,a,'mult',z)
                second=fit(X,squared_L,squared_U,[0]*3,'add',Q(0))
                check((first is None)==(second is None),'Squared feasibility correspondence')
                if first is not None:
                    expected=[Q(v['coefficient'])**2*z**v['sqrt_t_power'] for v in first['response']]
                    check(expected==[Q(v['coefficient']) for v in second['response']],'Squared prefix-max lift')
                squared_comparisons+=1
            # Interior midpoint tests and exact endpoint event removals.
        es=events(vals,vals,a,metric)
        for b,c in zip(es,es[1:]):
            z=(b+c)/2
            check((fit(X,vals,vals,a,metric,z) is not None)==order_oracle(vals,vals,a,metric,z),'Interior disagreement')
            comparisons+=1
    # Exact symbolic non-square boundary: t=2, gamma=sqrt(2).
    w=cases[0]['multiplicative_boundary_witness']
    check(any(y['sqrt_t_power']==1 for y in w['response']),'Non-square response representation')
    # Boundary tie is required: the latent peak cannot become strictly increasing.
    check(cases[0]['additive_boundary_witness']['response'][1]==cases[0]['additive_boundary_witness']['response'][2],'Expected plateau')
    # Bad negative original intervals must not be silently clipped at zero.
    try:events([Q(-1)],[Q(0)],[1],'add')
    except AssertionError:pass
    else:raise AssertionError('Negative input accepted')
    receipt={'status':'PASS_SYNTHETIC_ONLY','known_fixtures':cases,'weak_order_comparisons':comparisons,'squared_interval_equivalence_checks':squared_comparisons,'first_witness_vs_full_system':{'first_duplicate_order_epsilon':'1/2','full_epsilon':'3/2','first_duplicate_order_t':'2','full_t':'4'},'permanent_fixed_orders_checked':True,'unknown_upper_and_zero_boundary_checked':True,'non_square_t_symbolic_response_checked':True,'no_observed_source_loaded':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out=Path(__file__).with_name('synthetic_receipt_v2.json')
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':receipt['status'],'weak_order_comparisons':comparisons,'source_sha256':receipt['source_sha256']}))


if __name__=='__main__':main()
