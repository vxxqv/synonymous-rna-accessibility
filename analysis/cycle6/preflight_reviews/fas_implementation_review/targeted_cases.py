"""Additional synthetic full-design and redundant-feature adversarial checks."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
from adversarial_check import check_certificate, require

HERE=Path(__file__).resolve().parent
target=HERE/'candidate_v1/ordinal_geometry.py'
spec=importlib.util.spec_from_file_location('target',target)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
out=HERE/'targeted_results.json';require(not out.exists(),'No overwrite')
records=[]
states=list(itertools.product(range(3),range(3),range(2),range(2),range(2)))
design=m.features(states)
require(len(design)==72 and len(design[0])==8,'Design shape')
lower=[F(0)]*72;upper=[None]*72
for state,value in [((0,0,0,0,0),2),((0,1,0,0,0),0),((1,0,0,0,0),0),((1,1,0,0,0),2)]:
 i=states.index(state);lower[i]=F(value);upper[i]=F(value)
for perm in [list(range(72)),list(reversed(range(72))),list(range(1,72,2))+list(range(0,72,2))]:
 x=[design[i] for i in perm];l=[lower[i] for i in perm];u=[upper[i] for i in perm]
 r=m.order_certificate(x,l,u);require(r['status']=='infeasible','Hidden unknowns erased exact XOR obstruction')
 check_certificate(x,l,u,r);records.append({'case':'full72_exact_negative_with68unknown','dual':r['dual']})
q2=list(itertools.product(range(2),repeat=2))
for x in [[[1,a,b,a+b,0,2*a] for a,b in q2],[[1,1,0,0] for _ in q2]]:
 for vals,expected in [([0,1,2,3],'feasible' if len(set(map(tuple,x)))>1 else 'infeasible'),
                       ([2,0,0,2],'infeasible'),([1,1,1,1],'feasible')]:
  r=m.order_certificate(x,vals,vals);require(r['status']==expected,'Redundant-feature status mismatch')
  check_certificate(x,vals,vals,r);records.append({'case':'redundant_and_duplicate_columns','status':r['status']})
# Extremely small/large rational intervals do not become floating response thresholds.
for scale in [F(1,10**80),F(10**80),F(1)]:
 for values,want in [([0,1,2,3],'feasible'),([2,0,0,2],'infeasible')]:
  v=[scale*a for a in values];x=[[1,a,b] for a,b in q2]
  r=m.order_certificate(x,v,v);require(r['status']==want,'Rational response scale changed ordinal decision')
  check_certificate(x,v,v,r);records.append({'case':'exact_response_scale','scale':str(scale),'status':r['status']})
edges=[(0,1),(1,2)]
for e in [edges,list(reversed(edges)),edges+[(0,0),(0,1)],iter(edges)]:
 r=m.bounded_pairs([2,0,2],[2,None,2],e)
 require(r[1]['bottleneck_upper']=='2' and not r[1]['certified_nonempty'],'Unknown bridge or edge representation false exclusion')
 records.append({'case':'unknown_bridge_edge_representations','status':'no_exclusion'})
payload={'status':'PASS','case_count':len(records),'cases':records,
         'candidate_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
         'scope':'Synthetic rational values only, no biological inputs'}
out.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'PASS','cases':len(records),'sha256':hashlib.sha256(out.read_bytes()).hexdigest()}))
