"""Independent weak-order oracle for robustness search. Synthetic data only."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
import time

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
CORE=ROOT/'outputs/synonymous-rna-accessibility/analysis/cycle6/fas_pilot/ordinal_geometry.py'
EXPECTED='61c1f68f94bc57abf970b128f16f26aad44a1a9314eb0293f54dc288f4e5cb50'

def need(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def weak_orders():
    # Four central lines a=0,b=0,a=b,a=-b partition slope space.
    # Integer representatives in [-2,2]^2 hit all eight rays, eight sectors and origin.
    orders=set()
    for a,b in product(range(-2,3),repeat=2):
        latent=(0,a,b,a+b);levels=sorted(set(latent))
        orders.add(tuple(levels.index(v) for v in latent))
    need(len(orders)==17,'Weak-order arrangement coverage')
    return sorted(orders)

def order_radius(order,L,U,mask,metric):
    minimum=Q(0) if metric=='additive' else Q(1)
    # Greedy monotone interval feasibility requires every earlier lower <= later upper.
    for i,ri in enumerate(order):
        for j,rj in enumerate(order):
            if ri>rj or U[j] is None or L[i]<=U[j]:continue
            degree=mask[i]+mask[j]
            if degree==0 or metric!='additive' and U[j]==0:return math.inf
            if metric=='additive':candidate=(L[i]-U[j])/degree
            elif degree==1:candidate=(L[i]/U[j])**2
            else:candidate=L[i]/U[j]
            minimum=max(minimum,candidate)
    return minimum

def direct_oracle(orders,L,U,mask,metric):
    return min(order_radius(order,L,U,mask,metric) for order in orders)

def certificate_check(cert,X,L,U):
    n=len(L);forced={(i,j) for i in range(n) for j in range(n) if U[i] is not None and U[i]<L[j]}
    need(cert['forced_orders']==len(forced) and cert['exact_checked'] is True,'Certificate metadata')
    if cert['status']=='feasible':
        beta=list(map(Q,cert['beta']));latent=[sum(v*b for v,b in zip(row,beta)) for row in X]
        response=list(map(Q,cert['response']))
        need(list(map(Q,cert['latent']))==latent,'Latent witness')
        need(all(L[i]<=response[i] and (U[i] is None or response[i]<=U[i]) for i in range(n)),'Response interval witness')
        need(all(response[i]<=response[j] for i in range(n) for j in range(n) if latent[i]<=latent[j]),'Weak monotone witness')
        need(all(latent[j]-latent[i]>=1 for i,j in forced),'Unit strict margins')
    else:
        need(cert['status']=='infeasible','Unknown status')
        weights=[];total=[Q(0)]*len(X[0])
        for d in cert['dual']:
            i,j=d['lower_vertex'],d['higher_vertex'];w=Q(d['weight'])
            need((i,j) in forced and w>0,'Invalid dual support')
            row=[b-a for a,b in zip(X[i],X[j])]
            need(row==d['difference'],'Dual difference mismatch')
            total=[a+w*b for a,b in zip(total,row)];weights.append(w)
        need(sum(weights)==1 and all(v==0 for v in total),'Dual cancellation')

def independent_fields(L,U,mask,t,metric):
    if metric=='additive':
        return [max(Q(0),l-mask[i]*t) for i,l in enumerate(L)], [None if u is None else u+mask[i]*t for i,u in enumerate(U)]
    return [l*l/(t if mask[i] else 1) for i,l in enumerate(L)], [None if u is None else u*u*(t if mask[i] else 1) for i,u in enumerate(U)]

def audit_one(author,core,X,orders,L,U,mask,metric):
    expected=direct_oracle(orders,L,U,mask,metric)
    result=author.radius(core,X,L,U,mask,metric)
    if expected==math.inf:
        need(result['status']=='infinite','Expected infinite radius')
        last=Q(result['breakpoints'][-1]);lo,hi=independent_fields(L,U,mask,last,metric)
        certificate_check(result['permanent_dual'],X,lo,hi)
        for d in result['permanent_dual']['dual']:
            i,j=d['lower_vertex'],d['higher_vertex']
            need(mask[i]+mask[j]==0 or metric!='additive' and U[i]==0,'Nonpermanent infinity support')
    else:
        need(result['status']=='finite' and Q(result['parameter_exact'])==expected,'Threshold mismatch')
        lo,hi=independent_fields(L,U,mask,expected,metric)
        need(result['boundary_lower']==list(map(str,lo)) and result['boundary_upper']==[None if u is None else str(u) for u in hi],'Boundary fields')
        certificate_check(result['feasible_at_boundary'],X,lo,hi)
        if metric!='additive':
            a,b=Q(result['gamma_enclosure']['lower']),Q(result['gamma_enclosure']['upper'])
            need(a*a<=expected<=b*b,'Exact square-root enclosure')
        minimum=Q(0) if metric=='additive' else Q(1)
        if expected>minimum:
            cert=result['sharp_lower_certificate'];need(cert is not None,'Missing lower dual')
            # Check supported comparisons at a point before the boundary, then their exact disappearance time.
            prior=Q(result['breakpoints'][result['breakpoints'].index(str(expected))-1])
            l,u=independent_fields(L,U,mask,prior,metric);certificate_check(cert,X,l,u)
            for d in cert['dual']:
                i,j=d['lower_vertex'],d['higher_vertex'];a=mask[i]+mask[j]
                if a==0 or metric!='additive' and U[i]==0:continue
                death=(L[j]-U[i])/a if metric=='additive' else (L[j]/U[i])**(2 if a==1 else 1)
                need(death>=expected,'Dual support fails below sharp boundary')
    return result['status']

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    need(not args.out.exists(),'Refusing existing result');started=time.perf_counter()
    source=ROOT/'work/research_cycle7/robustness_before_import_guard_fix.py';raw=source.read_bytes()
    need(sha(raw)==EXPECTED,'Author source changed before review')
    snapshot=BASE/'original_robustness.py';need(not snapshot.exists(),'Refusing snapshot overwrite');snapshot.write_bytes(raw)
    proto=ROOT/'work/research_cycle7/robustness_protocol.md'
    (BASE/'original_protocol.md').write_bytes(proto.read_bytes())
    author=load(snapshot,'reviewed_robustness');core=load(CORE,'verified_core')
    X=[[1,0,0],[1,1,0],[1,0,1],[1,1,1]];orders=weak_orders()
    counts={'finite':0,'infinite':0,'point_systems':0,'interval_systems':0}
    for vals in product((Q(0),Q(1,2),Q(1)),repeat=4):
        for mask in product((0,1),repeat=4):
            for metric in ('additive','multiplicative_squared'):
                status=audit_one(author,core,X,orders,list(vals),list(vals),mask,metric)
                counts[status]+=1;counts['point_systems']+=1
    interval_alphabet=[(Q(0),Q(0)),(Q(1,2),Q(1)),(Q(0),None)]
    for raw_intervals in product(interval_alphabet,repeat=4):
        L=[a for a,b in raw_intervals];U=[b for a,b in raw_intervals]
        for mask in ((0,1,1,1),(0,0,0,0),(1,1,1,1)):
            for metric in ('additive','multiplicative_squared'):
                status=audit_one(author,core,X,orders,L,U,mask,metric)
                counts[status]+=1;counts['interval_systems']+=1
    # Generate 72 synthetic rows only, with no sequence or count-file input.
    states=list(product(range(3),range(3),range(2),range(2),range(2)))
    data={'vertices_data':[{'codon_states':list(s),'lower':'1','upper':'1','scenario_ratios':['1']*27} for s in states]}
    inventory=author.analyze(core,data)
    need(len(inventory['all_single_vertex_deletions'])==72,'Row-deletion coverage')
    for i,row in enumerate(inventory['all_single_vertex_deletions']):
        need(row['omitted_vertex']==i and row['retained_vertex_ids']==[j for j in range(72) if j!=i],'Row-deletion indexing')
        need(row['certificate']['status']=='feasible','Constant deletion witness')
    prep=inventory['all_single_preparation_deletions'];need(len(prep)==12,'Preparation-deletion coverage')
    for index,row in enumerate(prep):
        kind='output' if index<9 else 'input';omit=index if index<9 else index-9
        expected=[3*b+t for b in range(9) for t in range(3) if (b if kind=='output' else t)!=omit]
        need(row['kind']==kind and row['omitted_zero_based']==omit and row['retained_scenarios']==expected,'Preparation indexing')
        need(row['results']['additive']['parameter_exact']=='0' and row['results']['multiplicative_squared']['parameter_exact']=='1','Deletion known radii')
    # Preserve a precise pre-freeze-execution guard counterexample with no real result file.
    marker=BASE/'executed_before_freeze.txt';poison=BASE/'synthetic_core_marker.py'
    poison.write_text('from pathlib import Path\nPath('+repr(str(marker))+').write_text("executed")\n',encoding='utf-8')
    process=subprocess.run([sys.executable,str(snapshot),'--core',str(poison),'--out',str(BASE/'must_not_exist.json'),
                            '--result',str(BASE/'nonexistent_result.json'),'--freeze',str(BASE/'nonexistent_freeze.json')],capture_output=True,text=True)
    need(process.returncode!=0 and marker.exists() and not (BASE/'must_not_exist.json').exists(),'Guard reproducer unexpectedly changed')
    unknown_data={'vertices_data':[dict(row) for row in data['vertices_data']]}
    unknown_data['vertices_data'][0]={'codon_states':list(states[0]),'lower':'0','upper':None,'scenario_ratios':[None]+['1']*26}
    try:author.analyze(core,unknown_data)
    except Exception as error:unknown_failure={'type':type(error).__name__,'message':str(error)}
    else:unknown_failure=None
    # Diagnostic public-API input: a mask that compares equal to binary integers but is not integral-typed.
    try:
        floatmask=author.radius(core,[[1],[1]],[Q(1),Q(3,5)],[Q(1),Q(3,5)],[0.0,1.0],'multiplicative_squared')
        floatmask_result={'status':floatmask['status'],'parameter':floatmask.get('parameter_exact')}
    except Exception as error:floatmask_result={'type':type(error).__name__,'message':str(error)}
    huge=author.sqrt_enclosure(Q(10**40+1))
    result={'status':'NUMERICS_PASS_WITH_GUARD_BLOCKER','source_sha256':EXPECTED,
            'protocol_sha256':sha(proto.read_bytes()),'core_sha256':sha(CORE.read_bytes()),
            'weak_orders':17,'counts':counts,'row_deletions':72,'preparation_deletions':12,
            'guard_blocker':{'unverified_core_executed_before_freeze_gate':marker.exists(),'returncode':process.returncode,
                             'stderr':process.stderr},'unknown_preparation_diagnostic':unknown_failure,
            'float_mask_diagnostic':floatmask_result,'large_gamma_display_diagnostic':huge,
            'real_outcomes_read':False,'runtime_seconds':time.perf_counter()-started,
            'verifier_sha256':sha(Path(__file__).read_bytes())}
    args.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','counts','runtime_seconds','real_outcomes_read')}))

if __name__=='__main__':main()
