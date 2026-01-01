"""Compare synthetic primary results with the frozen independent reference."""
import sys
if not __debug__:raise RuntimeError('Disabled assertions refused before input/output')
sys.dont_write_bytecode=True
import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess
import time
import reference as ref

HELD='7f901518cfc196de5d6d24c10f6f1b336255063ca8bd22d544accb6ef7b940f3'
LABELS={'joint_connected':'jointly_connected','joint_only':'joint_only_obstruction','individual_disconnected':'single_trait_disconnection','ineligible_endpoint':'ineligible_endpoints'}
GRAPH=[[1,2],[0,3],[0,3],[1,2]]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--primary',required=True);parser.add_argument('--out',required=True);args=parser.parse_args()
    primary=Path(args.primary);out=Path(args.out);here=Path(__file__).resolve().parent
    if out.exists():raise FileExistsError(out)
    before={p.name:sha(p) for p in [primary,here/'reference.py',here/'primary_probe.py',Path(__file__)]}
    if sha(primary)!=HELD:raise ValueError('Held source changed')
    cases=[];expected=[];kinds=[];start=time.perf_counter()
    def add(case,value,kind):cases.append(case);expected.append(value);kinds.append(kind)
    for flat in product((-1,0,1),repeat=8):
        values=[list(flat[2*i:2*i+2]) for i in range(4)]
        for endpoints in [(0,3),(1,2)]:
            add({'operation':'classify','graph':GRAPH,'values':values,'endpoints':endpoints},
                ref.classify(dict(enumerate(values)),ref.SQUARE_EDGES,*endpoints),'point')
    options=[None,[-1,-1],[0,0],[1,1],[-1,1]]
    for flat in product(options,repeat=4):
        values=[[[0,0],[0,0]],[flat[0],flat[1]],[flat[2],flat[3]],[[0,0],[0,0]]]
        add({'operation':'interval','graph':GRAPH,'values':values,'endpoints':[0,3],'thresholds':[0,0]},
            ref.bounds(dict(enumerate(values)),ref.SQUARE_EDGES,0,3,(0,0)),'interval')
    for values in product((-1,0,1),repeat=4):
        add({'operation':'projection','values':values},ref.interaction(values),'projection')
    for target in range(4):
        ids=[v for v in range(4) if v!=target]
        for values in product((-1,0,1),repeat=3):
            known=dict(zip(ids,values));reference=ref.local_predictions({k:(v,v) for k,v in known.items()},target)
            add({'operation':'prediction','known':known,'target':target},reference,'prediction')
    for a,b,c,d in product((-1,0,1),repeat=4):
        values=[[1,1],[a,c],[b,d],[1,1]]
        add({'operation':'geometry','graph':GRAPH,'values':values,'endpoints':[0,3]},values,'geometry')
    replicates=[[1]*7,[0,1,2,0,2,1,0],[2,1,0,2,0,1,2],[1]*7]
    add({'operation':'preparations','graph':GRAPH,'values':replicates,'endpoints':[0,3]},
        ref.preparation_summary(dict(enumerate(replicates)))['preparation_statuses'],'preparations')
    child=subprocess.run([sys.executable,str(here/'primary_probe.py'),str(primary),HELD],
                         input=json.dumps(cases),capture_output=True,text=True,check=True)
    if child.stderr:raise RuntimeError('Unexpected primary stderr: '+child.stderr)
    observed=json.loads(child.stdout)
    if len(observed)!=len(cases):raise ValueError('Output truncation')
    counts={};boundary_probes=0
    for case,truth,result,kind in zip(cases,expected,observed,kinds):
        counts[kind]=counts.get(kind,0)+1
        if kind=='point':
            assert result['classification']==LABELS[truth['status']],(case,truth,result)
            if truth['status']!='ineligible_endpoint':
                assert result['individual_connectivity']==truth['individual'] and result['joint_connectivity']==truth['joint']
        elif kind=='interval':
            for key in ['individual_pessimistic','individual_optimistic','joint_pessimistic','joint_optimistic']:assert truth[key]==result[key],(case,key,truth,result)
            assert truth['joint_only_certificate']==(result['status']=='certified_joint_only_obstruction')
        elif kind=='projection':assert result['contrast']==truth['contrast'] and result['projection']==truth['projection']
        elif kind=='prediction':
            assert result['additive']==truth['additive'][0]
            assert result['adjacent_mean']==truth['neighbor_mean'][0]
            assert result['three_corner_mean']==truth['three_mean'][0]
        elif kind=='preparations':assert [x['classification'] for x in result]==[LABELS[x] for x in truth]
        elif kind=='geometry':
            for x,y in product((-3.,-1.,-.5,0.,.5,1.),repeat=2):
                indices=[]
                for coordinate,axis in zip((x,y),result['axes']):
                    matches=[]
                    for i,cell in enumerate(axis):
                        lo=None if cell['lower'] is None else float(Fraction(cell['lower']))
                        hi=float(Fraction(cell['upper']))
                        inside=coordinate==lo if cell['kind']=='atom' else (lo is None or coordinate>lo) and coordinate<hi
                        if inside:matches.append(i)
                    assert len(matches)==1,(coordinate,axis,matches)
                    indices.append(matches[0])
                matches=[cell for cell in result['cells'] if cell['axis_indices']==indices]
                assert len(matches)==1
                label=ref.classify(dict(enumerate(truth)),ref.SQUARE_EDGES,0,3,(x,y))['status']
                assert matches[0]['classification']==LABELS[label]
                boundary_probes+=1
    assert {p.name:sha(p) for p in [primary,here/'reference.py',here/'primary_probe.py',Path(__file__) ]}==before
    out.mkdir(parents=True)
    result={'status':'passed','counts':counts,'threshold_boundary_probes':boundary_probes,
      'total_primary_invocations':len(cases),'input_hashes_before_after_equal':before,
      'independence':'Reference imports no primary implementation. Primary API executes in an isolated probe process after reference freeze.',
      'scope':'Synthetic values only; no assay archive, arm values or observed result files opened.',
      'elapsed_seconds':time.perf_counter()-start}
    (out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
