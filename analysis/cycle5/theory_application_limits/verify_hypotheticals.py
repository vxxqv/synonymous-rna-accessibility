"""Exact checks of the report's hypothetical fields, with no data input."""
if not __debug__:raise SystemExit('Refuse disabled assertions')
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

def projection(y):
    h=(1,-1,-1,1);c=sum(sign*v for sign,v in zip(h,y))
    return c,[v-sign*c/4 for sign,v in zip(h,y)]
def effects(traits):
    return [((y[2]+y[3]-y[0]-y[1])/2,(y[1]+y[3]-y[0]-y[2])/2) for y in traits]
def determinant(rows):return rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0]
def hamming(a,b):return sum(x!=y for x,y in zip(a,b))

checks={}
e=F(3,7);y=[F(0),-2*e,-2*e,F(0)];c,p=projection(y)
assert c==4*e and p==[-e]*4 and max(abs(a-b) for a,b in zip(y,p))==e
assert min(max(-y[i],F(0)) for i in (1,2))==2*e
checks['constant_projection_sharp_2epsilon']=True
for d in (F(1),F(1,10**12),F(-1,10**12)):
    rows=effects([[F(0),F(-1),F(1),F(0)],[F(0),F(1),F(-1),d]])
    assert determinant(rows)==d
checks['arbitrarily_small_noise_rank_jump']=True
z=[F(0),F(2),F(1),F(3)];cubed=[x**3 for x in z]
c,p=projection(cubed);rows=effects([z,cubed])
assert c==18 and rows==[(F(1),F(2)),(F(10),F(17))] and determinant(rows)==-3
for s in range(4):
    for t in range(4):
        assert [v>=min(z[s],z[t]) for v in z]==[v>=min(cubed[s],cubed[t]) for v in cubed]
checks['strict_monotone_preserves_acceptance_changes_additivity_and_rank']=True
sample=['0000','0010','0011','1011','1111'];safe=['0000','1000','1100','1110','1111']
edges={(i,j) for i,a in enumerate(sample) for j,b in enumerate(sample) if i<j and hamming(a,b)==1}
assert edges=={(0,1),(1,2),(2,3),(3,4)}
def f(x):return sum(map(int,x[:2]))-sum(map(int,x[2:]))
assert [f(x) for x in sample]==[0,-1,-2,-1,0]
assert all(hamming(a,b)==1 for a,b in zip(safe,safe[1:])) and min(f(x) for x in safe)>=0
checks['sampled_path_breaks_product_bound_full_cube_has_safe_route']=True
root=Path(__file__).resolve().parent
receipt={'status':'passed','checks':checks,'arithmetic':'exact fractions and integer Hamming graphs',
    'scientific_data_files_read':0,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'python':sys.version,'scope':'Displayed hypothetical examples only, not empirical application or a new proof of C5'}
with (root/'synthetic_checks.json').open('x',encoding='utf-8') as stream:json.dump(receipt,stream,indent=2);stream.write('\n')
print(json.dumps({'status':'passed','checks':len(checks)}))
