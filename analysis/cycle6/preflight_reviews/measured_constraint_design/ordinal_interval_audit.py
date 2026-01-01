"""Exact synthetic Fourier-Motzkin checks of a finite interval/ordinal contract.

No assay data, production solver or scientific model is imported.
"""
if not __debug__:raise RuntimeError('Optimization refused')
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def solve(rows,d):
 """Rows (a,r) mean a.x>=r. Return a rational witness or None."""
 if d==0:return () if all(r<=0 for a,r in rows) else None
 pos=[(a,r) for a,r in rows if a[-1]>0]
 neg=[(a,r) for a,r in rows if a[-1]<0]
 reduced=[(a[:-1],r) for a,r in rows if a[-1]==0]
 for a,r in pos:
  for b,s in neg:
   reduced.append((tuple(-b[-1]*a[j]+a[-1]*b[j] for j in range(d-1)),-b[-1]*r+a[-1]*s))
 prefix=solve(reduced,d-1)
 if prefix is None:return None
 lower=[(r-sum(a[j]*prefix[j] for j in range(d-1)))/a[-1] for a,r in pos]
 upper=[(r-sum(a[j]*prefix[j] for j in range(d-1)))/a[-1] for a,r in neg]
 lo=max(lower) if lower else None;hi=min(upper) if upper else None
 require(lo is None or hi is None or lo<=hi,'Elimination recovery contradiction')
 z=Q(0)
 if lo is not None:z=max(z,lo)
 if hi is not None:z=min(z,hi)
 answer=prefix+(z,)
 require(all(sum(c*v for c,v in zip(a,answer))>=r for a,r in rows),'Invalid rational witness')
 return answer
def inequalities(X,intervals):
 require(len(X)==len(intervals)>0 and all(len(x)==len(X[0]) for x in X),'Invalid feature dimensions')
 require(all(U is None or L<=U for L,U in intervals),'Empty input interval')
 # None is an explicit +infinity sentinel for upper bounds only.
 pairs=[(i,j) for i in range(len(X)) for j in range(len(X)) if intervals[i][1] is not None and intervals[i][1]<intervals[j][0]]
 return pairs,[(tuple(Q(y-x) for x,y in zip(X[i],X[j])),Q(1)) for i,j in pairs]
def fitted_values(X,intervals,b):
 scores=[sum(Q(x)*v for x,v in zip(row,b)) for row in X]
 values=[max(intervals[k][0] for k in range(len(X)) if scores[k]<=scores[i]) for i in range(len(X))]
 require(all(L<=y and (U is None or y<=U) for y,(L,U) in zip(values,intervals)),'Interval lifting failed')
 require(all(values[i]<=values[j] for i in range(len(X)) for j in range(len(X)) if scores[i]<=scores[j]),'Monotone lifting failed')
 return scores,values
def weak_orders(n):
 return [r for r in product(range(n),repeat=n) if set(r)==set(range(max(r)+1))]
def order_feasible(X,ranks):
 rows=[]
 for i in range(len(X)):
  for j in range(i):
   diff=tuple(Q(x-y) for x,y in zip(X[i],X[j]))
   if ranks[i]==ranks[j]:rows += [(diff,Q(0)),(tuple(-x for x in diff),Q(0))]
   elif ranks[i]>ranks[j]:rows.append((diff,Q(1)))
   else:rows.append((tuple(-x for x in diff),Q(1)))
 return solve(rows,len(X[0])) is not None
def interval_order_compatible(intervals,ranks):
 # Independent oracle: all prefix lower bounds must fit every upper bound
 # of the current tied group. No forced-pair construction is used here.
 for level in sorted(set(ranks)):
  prefix=max(intervals[k][0] for k,r in enumerate(ranks) if r<=level)
  ceilings=[intervals[k][1] for k,r in enumerate(ranks) if r==level and intervals[k][1] is not None]
  if ceilings and prefix>min(ceilings):return False
 return True
def main():
 matrices=[[(0,0),(1,0),(2,0)],[(0,0),(1,0),(0,1)],[(0,0),(0,0),(1,0)],[(0,0),(1,1),(-1,-1)]]
 grid=[(Q(l),Q(u)) for l in range(3) for u in range(l,3)]+[(Q(l),None) for l in range(3)]
 counts={'instances':0,'feasible':0,'infeasible':0,'dual_certificates':0,'weak_orders_per_matrix':13}
 for X in matrices:
  orders=[r for r in weak_orders(3) if order_feasible(X,r)]
  for intervals in product(grid,repeat=3):
   pairs,rows=inequalities(X,intervals);b=solve(rows,2)
   independent=any(interval_order_compatible(intervals,r) for r in orders)
   require((b is not None)==independent,'Forced-order LP disagrees with weak-order oracle')
   counts['instances']+=1
   if b is not None:
    counts['feasible']+=1;fitted_values(X,intervals,b)
   else:
    counts['infeasible']+=1;m=len(rows)
    dual=[]
    for j in range(2):
     col=tuple(a[j] for a,r in rows);dual += [(col,Q(0)),(tuple(-v for v in col),Q(0))]
    dual += [(tuple(Q(1) for _ in range(m)),Q(1)),(tuple(Q(-1) for _ in range(m)),Q(-1))]
    dual += [(tuple(Q(int(i==j)) for i in range(m)),Q(0)) for j in range(m)]
    lam=solve(dual,m)
    require(lam is not None and sum(lam)==1 and all(v>=0 for v in lam),'Missing dual certificate')
    require(all(sum(lam[i]*rows[i][0][j] for i in range(m))==0 for j in range(2)),'Dual cancellation fails')
    counts['dual_certificates']+=1
 require(counts['instances']==2916,'Wrong finite universe')
 # Minimal positive reciprocal sign pattern: an exact cancellation certificate.
 X=[(0,0),(1,0),(0,1),(1,1)]
 intervals=[(Q(1),Q(1)),(Q(0),Q(0)),(Q(0),Q(0)),(Q(1),Q(1))]
 pairs,rows=inequalities(X,intervals)
 require(solve(rows,2) is None,'Reciprocal-order pattern incorrectly additive')
 i=pairs.index((1,0));j=pairs.index((2,3));lam=[Q(0)]*len(rows);lam[i]=lam[j]=Q(1,2)
 require(all(sum(lam[k]*rows[k][0][a] for k in range(len(rows)))==0 for a in range(2)),'Q2 dual fails')
 # A wider Q2 envelope can remove the exact point contradiction.
 p,r=inequalities(X,[(Q(0),Q(1))]*4);b=solve(r,2)
 require(b is not None,'Overlapping Q2 intervals should allow a constant model');fitted_values(X,[(Q(0),Q(1))]*4,b)
 # Closed boundary and nonlinear plateau examples.
 plateau_X=[(0,),(1,),(2,)]
 plateau_I=[(Q(0),Q(0)),(Q(0),Q(0)),(Q(1),Q(1))]
 p,r=inequalities(plateau_X,plateau_I);b=solve(r,1)
 require(b is not None,'Nondecreasing plateau incorrectly rejected');fitted_values(plateau_X,plateau_I,b)
 touching=[(Q(0),Q(1)),(Q(1),Q(2)),(Q(1),Q(1))]
 p,r=inequalities(plateau_X,touching);require(not p,'Touching intervals forced an order')
 scores,values=fitted_values(plateau_X,touching,(Q(0),));require(values==[Q(1)]*3,'Constant/touching interval failure')
 zero_X=[(),()]
 p,r=inequalities(zero_X,[(Q(0),None),(Q(1),Q(2))]);b=solve(r,0)
 require(b==(),'Zero-dimensional constant model failed');fitted_values(zero_X,[(Q(0),None),(Q(1),Q(2))],b)
 p,r=inequalities(zero_X,[(Q(0),Q(0)),(Q(1),None)]);require(solve(r,0) is None,'Zero-dimensional disjoint intervals passed')
 try:inequalities([(0,)],[(Q(2),Q(1))])
 except RuntimeError:pass
 else:raise RuntimeError('Empty interval accepted')
 out={'status':'EXACT_INTERVAL_ORDINAL_SYNTHETIC_CHECKS_PASS','counts':counts,'Q2_forced_pairs':pairs,'Q2_Farkas_weights':[str(v) for v in lam],'plateau_requires_nondecreasing_not_strictly_increasing':True,'touching_intervals_and_constant_response_checked':True,'duplicate_features_and_latent_ties_checked':True,'unbounded_upper_intervals_checked':True,'zero_dimensional_and_empty_interval_guards_checked':True,'scope':'Synthetic mathematical feasibility only, not FAS observed data, confidence calibration or production solver validation','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 p=Path(__file__).with_name('ordinal_interval_checks_v3.json')
 if p.exists():raise RuntimeError('No overwrite')
 p.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(out,indent=2))

if __name__=='__main__':main()
