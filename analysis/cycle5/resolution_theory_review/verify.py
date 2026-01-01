"""Exact, finite, synthetic resolution checks; no external-data input option."""
if not __debug__:
    raise SystemExit('Disabled assertions refused before any input/output')
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import sys
import time
sys.dont_write_bytecode=True

EDGES=((0,1),(0,2),(1,3),(2,3))

def connected(allowed,s,t):
    if s not in allowed or t not in allowed:return False
    reached={s};frontier=[s]
    while frontier:
        u=frontier.pop()
        for a,b in EDGES:
            if u not in (a,b):continue
            v=b if u==a else a
            if v in allowed and v not in reached:reached.add(v);frontier.append(v)
    return t in reached

def classify(field,endpoints=(0,3),threshold=None):
    s,t=endpoints
    tau=threshold if threshold is not None else tuple(min(field[s][k],field[t][k]) for k in (0,1))
    accepted=[{v for v in range(4) if field[v][k]>=tau[k]} for k in (0,1)]
    individual=[connected(a,s,t) for a in accepted]
    joint=connected(accepted[0]&accepted[1],s,t)
    return 'joint_connected' if joint else 'joint_only' if all(individual) else 'individual_disconnected'

def box_certificate(field,eps):
    lower=[[x-eps for x in row] for row in field]
    upper=[[x+eps for x in row] for row in field]
    tau=[min(lower[0][k],lower[3][k]) for k in (0,1)]
    pess=[{v for v in range(4) if lower[v][k]>=tau[k]} for k in (0,1)]
    optimistic={v for v in range(4) if all(upper[v][k]>=tau[k] for k in (0,1))}
    return all(connected(a,0,3) for a in pess) and not connected(optimistic,0,3)

def grid(field,eps):
    choices=sorted({-eps,F(0),eps});counts=Counter()
    for errors in product(choices,repeat=8):
        perturbed=[[field[v][k]+errors[2*v+k] for k in (0,1)] for v in range(4)]
        counts[classify(perturbed)]+=1
    return dict(sorted(counts.items()))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    if args.out.exists():raise FileExistsError('Refuse existing output')
    source=Path(__file__);before=hashlib.sha256(source.read_bytes()).hexdigest();start=time.perf_counter()
    cases=[]
    for m,eps in [(F(1,2),F(1)),(F(1),F(1)),(F(2),F(1)),(F(3),F(1)),(F(1),F(1,2)),(F(3,2),F(1,2)),(F(1),F(0))]:
        field=[[F(0),F(0)],[m,-m],[-m,m],[F(0),F(0)]]
        counts=grid(field,eps);predicted=m>2*eps
        assert (set(counts)=={'joint_only'})==predicted
        assert box_certificate(field,eps)==predicted
        if not predicted:
            witness=[[-eps,-eps],[m-eps,-m+eps],[-m,m],[-eps,-eps]]
            assert classify(witness)=='joint_connected'
            assert all(abs(witness[v][k]-field[v][k])<=eps for v in range(4) for k in (0,1))
        assert classify(field,(1,2))=='joint_connected'
        for c1,c2 in product((F(-5),F(0),F(5)),repeat=2):
            shifted=[[row[0]+c1,row[1]+c2] for row in field]
            assert classify(shifted)==classify(field)
        cases.append({'m':str(m),'epsilon':str(eps),'universal_prediction':predicted,'grid_counts':counts,'fixed_low_range_certificate':box_certificate(field,eps)})
    # Positive equality is sufficient under weak acceptance, despite no strict margin.
    weak=[[F(0),F(0)],[F(2),F(-3)],[F(-3),F(2)],[F(0),F(0)]]
    weak_counts=grid(weak,F(1));assert set(weak_counts)=={'joint_only'}
    # Fixed-low certificate and recomputed endpoint-relative quantifier differ.
    asymmetric=[[F(0),F(0)],[F(1),F(-3)],[F(-3),F(1)],[F(0),F(0)]]
    witness=[[F(1),F(1)],[F(0),F(-4)],[F(-4),F(0)],[F(1),F(1)]]
    assert box_certificate(asymmetric,F(1))
    assert classify(witness)=='individual_disconnected'
    assert all(abs(witness[v][k]-asymmetric[v][k])<=1 for v in range(4) for k in (0,1))
    # The comparison bound is attained in each direction at a nonendpoint.
    for sign in (-1,1):
        baseline=F(7,3);eps=F(2,5)
        original_gap=baseline-min(F(0),F(1))
        changed=baseline+sign*eps-min(-sign*eps,F(1)-sign*eps)
        assert changed-original_gap==sign*2*eps
    after=hashlib.sha256(source.read_bytes()).hexdigest();assert before==after
    result={'status':'passed','assertion_tested':'Exact rational ternary-grid checks and declared witnesses, not exhaustive real-box verification',
        'symmetric_cases':cases,'weak_positive_boundary_grid_counts':weak_counts,
        'fixed_low_range_counterexample':{'field':[[str(x) for x in row] for row in asymmetric],
        'perturbation_witness':[[str(x) for x in row] for row in witness],
        'range_certificate':True,'endpoint_relative_witness_status':classify(witness)},
        'total_ternary_grid_fields':sum(sum(c['grid_counts'].values()) for c in cases)+sum(weak_counts.values()),
        'source_sha256_before':before,'source_sha256_after':after,
        'argv':sys.argv,'python':sys.version,'arithmetic':'fractions.Fraction, exact comparisons, no tolerance',
        'elapsed_seconds':time.perf_counter()-start,'biological_files_read':0,
        'proof_dependency':'Universal continuum claim is proved in review.md independently of finite checks.'}
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/'receipt.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':'passed','grid_fields':result['total_ternary_grid_fields'],'source_sha256':before}))

if __name__=='__main__':main()
