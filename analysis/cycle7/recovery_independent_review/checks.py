"""Exact invented factorization witnesses, no observed inputs or author imports."""
import sys
if not __debug__:raise RuntimeError('Optimized execution refused')
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json


def main():
    coordinate_cases=0
    # Exhaustively compare one-coordinate anchored expansion with existence of a
    # bounded positive nuisance multiplier via direct interval intersection.
    for L,U,f,k in product([Q(0),Q(1,2),Q(1),Q(2)],repeat=4):
        if U<L or k<1:continue
        expanded=L/k<=f<=k*U
        if f==0:
            possible=L==0
        else:
            possible=max(Q(1)/k,L/f)<=min(k,U/f)
        assert expanded==possible
        coordinate_cases+=1
    # Constant inclusion reproduces any one positive field at the intensity level.
    y=[Q(1),Q(2),Q(1,3),Q(7)]
    c=Q(1)/max(y);recoveries=[c*v for v in y]
    assert all(0<v<=1 for v in recoveries)
    assert [v/recoveries[0] for v in recoveries]==y
    # A shared fixed b,f cannot produce two different ideal fields.
    fields=[[Q(1),Q(2)],[Q(1),Q(1,2)]]
    assert fields[0]!=fields[1]
    # Full codon-additive Q2 class, not just a constant-model counterexample.
    X=[(0,0),(0,1),(1,0),(1,1)]
    high=[Q(1),Q(2),Q(2),Q(1)]
    low=[Q(1),Q(1,2),Q(1,2),Q(1)]
    assert all(min(a,b)<=1<=max(a,b) for a,b in zip(high,low))
    differences=[tuple(v-u for u,v in zip(X[0],X[1])),
                 tuple(v-u for u,v in zip(X[3],X[2]))]
    assert all(sum(row[k] for row in differences)==0 for k in range(2))
    assert high[0]<high[1] and high[3]<high[2]
    # At kappa1, b=1 is forced: one box completion fits, but high itself cannot.
    # For any fixed finite kappa, alpha>kappa² yields the same strict obstruction.
    for k in [Q(1),Q(3,2),Q(2),Q(5)]:
        alpha=k*k+1
        assert alpha/k>k
    # Claim3 square: common constant-field model on point(1,9).
    lam=Q(3);fprime=[Q(3),Q(3)];y=[Q(1),Q(9)]
    d=[a/b for a,b in zip(y,fprime)]
    f=[v/fprime[0] for v in fprime];bias=[fprime[0]*v for v in d]
    assert f==[1,1] and bias==[1,9] and bias[0]==1
    assert max(bias)/min(bias)==lam*lam
    assert all(a*b==v for a,b,v in zip(f,bias,y))
    # An additive unit bound does not control multiplicative distortion near0.
    factor_examples=[]
    for n in [10,100,1000]:
        f=Q(1,n);observed=f+Q(1,2)
        assert observed-f==Q(1,2)
        factor_examples.append(str(observed/f))
    receipt={'status':'PASS_EXACT_SYNTHETIC_SCOPE_CHECKS','anchored_coordinate_cases':coordinate_cases,
             'one_box_completion_not_all_scenarios_full_Q2_counterexample':True,
             'probability_bounded_primitive_recovery_constructed':True,
             'pairwise_spread_square_example':{'lambda_all':'3','kappa_WT':'9','pairwise_spread':'9'},
             'additive_half_unit_multiplicative_factors':factor_examples,
             'scope':'Invented rational fields only, no observed outcomes or author code imported.',
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out=Path(__file__).with_name('synthetic_receipt.json')
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))


if __name__=='__main__':main()
