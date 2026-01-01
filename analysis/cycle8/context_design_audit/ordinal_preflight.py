"""Exact invented cases for the proposed interval-order model contract."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json

def constraints(X,intervals):
    out=[]
    for i,(li,ui) in enumerate(intervals):
        for j,(lj,uj) in enumerate(intervals):
            if uj is not None and li>uj:
                out.append({'higher':i,'lower':j,'a':[F(x)-F(y) for x,y in zip(X[i],X[j])]})
    return out

def completion(X,beta,intervals):
    z=[sum(F(x)*b for x,b in zip(row,beta)) for row in X]
    for c in constraints(X,intervals):
        if sum(a*b for a,b in zip(c['a'],beta))<1:raise ValueError('Invalid margin witness')
    positives=[v for p in intervals for v in p if v is not None and v>0]
    h=min(positives)/2 if positives else F(1);values=[None]*len(z)
    for level in sorted(set(z)):
        inds=[i for i,v in enumerate(z) if v==level]
        h=max(h,max(intervals[i][0] for i in inds))
        if any(intervals[i][1] is not None and h>intervals[i][1] for i in inds):raise ValueError('No monotone completion')
        for i in inds:values[i]=h
    return z,values

def dual(cs,weights):
    if len(cs)!=len(weights) or any(w<0 for w in weights) or sum(weights)!=1:raise ValueError('Invalid dual weights')
    if any(sum(w*c['a'][j] for w,c in zip(weights,cs)) for j in range(len(cs[0]['a']))):raise ValueError('Invalid dual cancellation')
    return True

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    if args.out.exists():raise FileExistsError(args.out)
    X=[[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
    # Checkerboard: both positive-sign corners exceed both negative-sign corners.
    ints=[(F(2),F(2)),(F(1),F(1)),(F(1),F(1)),(F(2),F(2))]
    cs=constraints(X,ints);weights=[F(1,2) if (c['higher'],c['lower']) in [(0,1),(3,2)] else F(0) for c in cs]
    dual(cs,weights)
    I=[[int(i==j) for j in range(4)] for i in range(4)]
    zi,yi=completion(I,[F(2),F(1),F(1),F(2)],ints)
    # The arbitrary common-order class passes one context even though additive fails.
    order=[1,2,0,3]
    if [ints[i][0] for i in order]!=sorted(p[0] for p in ints):raise ValueError('Common order example')
    # Touching closed intervals do not force a strict order; unknown has open zero.
    touch=[(F(1),F(2)),(F(2),F(3)),(F(0),None)]
    if constraints([[1],[1],[1]],touch):raise ValueError('False forced tie order')
    zt,yt=completion([[1],[1],[1]],[F(0)],touch)
    if yt!=[F(2)]*3:raise ValueError('Tie completion')
    # Every rational strictly feasible finite system admits global rescaling to margin one.
    small=[F(0),F(1,100),F(2,100),F(3,100)]
    first=[(x,x) for x in small]
    A=constraints(I,first);gaps=[sum(a*b for a,b in zip(c['a'],small)) for c in A]
    scaled=[v/min(gaps) for v in small]
    if any(sum(a*b for a,b in zip(c['a'],scaled))<1 for c in A):raise ValueError('Rescaling')
    # Two cell contexts reverse an order. Shared scores fail; independent slopes pass.
    two=[[1,0],[1,1]];a=[(F(1),F(1)),(F(2),F(2))];b=list(reversed(a))
    both=constraints(two,a)+constraints(two,b);dual(both,[F(1,2)]*2)
    za,ya=completion(two,[F(0),F(1)],a);zb,yb=completion(two,[F(1),F(-1)],b)
    # Widening never introduces a new strict pair requirement.
    old={(c['higher'],c['lower']) for c in cs};widened=[]
    for g in [F(1),F(5,4),F(3,2),F(2)]:
        w=[(l/g,None if u is None else u*g) for l,u in ints]
        new={(c['higher'],c['lower']) for c in constraints(X,w)}
        if not new<=old:raise ValueError('Widening strengthened constraints')
        widened.append({'gamma':str(g),'forced_pairs':len(new)})
    result={'status':'PASS','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'cases':['additive checkerboard exact normalized dual','saturated pair-model exact feasible witness',
        'unrestricted order not equivalent to additivity','touching intervals and latent ties','unknown open-zero positive completion',
        'strict finite gaps rescale to unit margin','shared versus cell-specific separation','all four fixed widening factors'],
      'checkerboard_constraints':cs,'checkerboard_dual':weights,'pair_model_values':yi,
      'tied_values':yt,'widening':widened,'biological_outcomes_read':False}
    args.out.mkdir(parents=True);(args.out/'receipt.json').write_text(json.dumps(result,indent=2,default=str)+'\n')
    print('PASS: 8 exact synthetic model-contract cases; no biological outcomes read.')
if __name__=='__main__':main()
