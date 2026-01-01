"""Synthetic exact-arithmetic claim checks, no source outcomes or author imports."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,hashlib

if not __debug__:raise RuntimeError('Optimized Python prohibited')
p=Path('work/research_cycle8/comt_independent_review');p.mkdir(parents=True,exist_ok=True)
dest=p/'shared_reference_checks.json'
if dest.exists():raise FileExistsError(dest)

def mean(x):return sum(x,F())/len(x)
def cov(x,y):return mean([(a-mean(x))*(b-mean(y)) for a,b in zip(x,y)])

rows=[{'R':F(r-d),'T':F(d-f),'L':F(h-r)} for r,d,h,f in product([-1,1],repeat=4)]
columns={k:[r[k] for r in rows] for k in ['R','T','L']}
covs={a+b:str(cov(columns[a],columns[b])) for a,b in [('R','T'),('R','L'),('T','L')]}
assert covs=={'RT':'-1','RL':'-1','TL':'0'}

# True biological RNA and heavy abundance both equal u, so ideal heavy/RNA load is constant.
# Only RNA measurement is perturbed by e. The low-fluorescence proxy is perfectly -u.
# Thus observed R=u+e,T=u,L=-e and a=1,b=1 predicts perfectly without variable true load.
toy=[]
for u,e in product([-1,1],[-2,2]):
    r,d,h,f=F(u+e),F(0),F(u),F(-u)
    toy.append({'R':r-d,'T':d-f,'L':h-r,'true_load':F(0)})
toy.sort(key=lambda x:x['R'])
assert [r['R'] for r in toy]==list(map(F,[-3,-1,1,3]))
assert [r['T'] for r in toy]==list(map(F,[-1,1,-1,1]))
isotonic=list(map(F,[-1,0,0,1]))
constant_mae=mean([abs(r['T']) for r in toy])
pava_mae=mean([abs(r['T']-y) for r,y in zip(toy,isotonic)])
joint_mae=mean([abs(r['T']-(r['R']+r['L'])) for r in toy])
assert (constant_mae,pava_mae,joint_mae)==(F(1),F(1,2),F(0))
out={'scope':'Synthetic algebra/claim counterexamples only; no biological values or author implementation imported',
     'independent_centered_raw_error_covariances':covs,
     'constant_true_load_positive_control':{'rows':[{k:str(v) for k,v in r.items()} for r in toy],
        'constrained_joint_coefficients':{'RNA_slope':'1','load_slope':'1','intercept':'0'},
        'constant_MAE':str(constant_mae),'nondecreasing_PAVA_MAE':str(pava_mae),'joint_MAE':str(joint_mae),
        'interpretation':'A capable joint model can improve solely through RNA measurement-error cancellation when true load is constant. No mechanism conclusion follows.'},
     'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
dest.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS_EXACT_SYNTHETIC_CHECKS','receipt_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()}))
