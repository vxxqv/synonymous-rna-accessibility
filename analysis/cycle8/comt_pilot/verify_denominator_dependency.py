"""Exact synthetic demonstration of why target-ID ascertainment cannot recover WT."""
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
P=Path(__file__).resolve().parent
def normalized(depth,retained_counts,target):
    if any(type(x)is not int or x<0 for x in [depth,*retained_counts]):raise ValueError('count domain')
    wt=depth-sum(retained_counts)
    if wt<0:raise ValueError('negative source WT')
    return F(2*target+1,2*wt+1)
def run():
    # Same accepted target and depth; only a background allele's RF retention changes.
    with_background=normalized(100,[10,30],10)
    without_background=normalized(100,[10],10)
    assert with_background==F(21,121) and without_background==F(21,181)
    assert with_background!=without_background
    # A shared WT is exactly the same random object, not independent per contrast.
    cancellation=[]
    for wt in (0,10,100):
        a=F(21,2*wt+1);b=F(41,2*wt+1)
        cancellation.append(b/a)
    assert cancellation==[F(41,21)]*3
    # Replacing it by separate references destroys that algebraic cancellation.
    assert F(41,201)/F(21,21)!=F(41,21)
    # Positive raw count does not guarantee positive background-subtracted representation.
    nc=F(21,121)
    assert with_background-nc==0
    assert without_background-nc<0
    # Source ROI2 four-sample filter convention is not prose strict majority.
    assert (2>=2) is True and (2>4/2) is False
    out={'status':'PASS','exact_checks':7,'target_ratio_background_retained':str(with_background),'target_ratio_background_rejected':str(without_background),'same_target_ID_but_different_normalized_value':True,'source_counts_read':False,'scope':'Synthetic dependency counterexample, shared-reference cancellation, nonpositive subtraction and exact-half filtering. Not an observed outcome.'}
    (P/'synthetic_dependency_receipt.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out))
if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions required')
    run()
