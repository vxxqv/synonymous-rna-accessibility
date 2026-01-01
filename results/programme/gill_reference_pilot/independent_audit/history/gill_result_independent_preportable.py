"""Independent saved-vector audit of the frozen Gill pilot, no workbook access.

Uses Decimal logarithms, conditional-mean factorial decomposition, the min-max
isotonic characterization, a Laplacian eigendecomposition, and threshold BFS.
No import from the fitted runner and no fitting choices or endpoint changes.
"""
from pathlib import Path
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations
from collections import deque
from datetime import datetime,timezone
import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/programme/gill_reference_pilot'
PROTOCOL='analysis/programme/retrospective/gill_reference_protocol_v1.json'
RUNNER='analysis/programme/retrospective/gill_reference_pilot.py'
MODELS=['constant','gc','additive','isotonic_additive','pairwise_half','heat','training_mean']
CS=['0.25','0.5','1.0']
TOL=2e-12

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def check_close(got,expected):
    got=np.asarray(got,dtype=float);expected=np.asarray(expected,dtype=float)
    if got.shape!=expected.shape or not np.isfinite(got).all() or not np.isfinite(expected).all():
        raise AssertionError('Nonfinite or shape-mismatched comparison')
    difference=float(np.max(abs(got-expected))) if got.size else 0.
    if difference>TOL:raise AssertionError(f'Numerical discrepancy {difference}')
    return difference

def conditional_effect(values,vertices,order):
    output=[]
    for vertex in vertices:
        value=Decimal(0)
        for selected in combinations(range(3),order):
            for size in range(order+1):
                for subset in combinations(selected,size):
                    matching=[y for x,y in zip(vertices,values) if all(x[k]==vertex[k] for k in subset)]
                    value+=(-1)**(order-size)*sum(matching)/len(matching)
        output.append(value)
    return output

def anchor(values,reference):return [v-values[reference] for v in values]

def score(counts,rna,dna,c):
    logtwo=Decimal(2).ln();pseudo=Decimal(c)
    return [sum((Decimal(counts[col][i])+pseudo).ln()/logtwo for col in rna)/len(rna)
            -sum((Decimal(counts[col][i])+pseudo).ln()/logtwo for col in dna)/len(dna)
            for i in range(32)]

def independent_isotonic(inputs,outputs):
    knots=sorted(set(inputs));weights=[inputs.count(k) for k in knots]
    totals=[sum(y for x,y in zip(inputs,outputs) if x==k) for k in knots]
    fit=[max(min(sum(totals[a:b+1])/sum(weights[a:b+1]) for b in range(i,len(knots)))
             for a in range(i+1)) for i in range(len(knots))]
    return [fit[knots.index(x)] for x in inputs]

def fit_independent(z,keys,vertices,reference,heat):
    parts=[conditional_effect(z,vertices,o) for o in range(4)]
    additive=[a+b for a,b in zip(parts[0],parts[1])]
    gc=[sum(b in 'GC' for b in key) for key in keys]
    meangc=sum(map(Decimal,gc))/32;mean=sum(z)/32
    slope=sum((Decimal(g)-meangc)*(v-mean) for g,v in zip(gc,z))/sum((Decimal(g)-meangc)**2 for g in gc)
    raw={'constant':[mean]*32,'gc':[mean+slope*(Decimal(g)-meangc) for g in gc],
         'additive':additive,'isotonic_additive':independent_isotonic(additive,z),
         'pairwise_half':[a+b/2 for a,b in zip(additive,parts[2])],
         'training_mean':z}
    result={name:list(map(float,anchor(value,reference))) for name,value in raw.items()}
    h=heat@np.array(list(map(float,z)));result['heat']=(h-h[reference]).tolist()
    return result

def independent_metrics(prediction,outcome,reference):
    error=np.array(prediction)-np.array(outcome)
    assert error[reference]==0
    mean=float(np.mean(error))
    return {'errors':error.tolist(),'rmse32':float(np.sqrt(np.dot(error,error)/32)),
            'rmse31':float(np.sqrt(np.dot(error,error)/31)),
            'mean_signed_error':mean,'centered_rmse':float(np.linalg.norm(error-mean)/math.sqrt(32))}

def independent_geometry(field,vertices,graph):
    capacity=np.full((32,32),np.nan);thresholds=[]
    for t in sorted(set(field),reverse=True):
        remaining={i for i,v in enumerate(field) if v>=t};n=len(remaining);sizes=[]
        while remaining:
            u=min(remaining);remaining.remove(u);component={u};todo=deque([u])
            while todo:
                for v in graph[todo.popleft()]:
                    if v in remaining:remaining.remove(v);component.add(v);todo.append(v)
            for u in component:
                for v in component:
                    if np.isnan(capacity[u,v]):capacity[u,v]=t
            sizes.append(len(component))
        thresholds.append({'threshold':t,'retained':n,'components':len(sizes),'component_sizes':sorted(sizes,reverse=True)})
    barriers={(i,j):min(field[i],field[j])-capacity[i,j] for i in range(32) for j in range(i+1,32)}
    d=[Decimal.from_float(float(v)) for v in field]
    energies={str(o):float(sum(v*v for v in conditional_effect(d,vertices,o))/32) for o in range(4)}
    return capacity,barriers,thresholds,energies

def main():
    if not __debug__:raise RuntimeError('Audit assertions require normal Python execution')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,default=BASE/'independent_audit')
    args=parser.parse_args()
    started=time.perf_counter();started_at=datetime.now(timezone.utc).isoformat();out=args.outdir.resolve()
    if not out.is_relative_to((ROOT/'results').resolve()):raise ValueError('Audit outputs must be inside repository results')
    if out.exists():raise FileExistsError('Use a fresh independent audit directory')
    rels=[PROTOCOL,RUNNER,'results/programme/retrospective/gill_amplicon_pilot_v2.json',
          'results/programme/gill_reference_pilot/prevalidation_gate.json']
    rels += ['results/programme/gill_reference_pilot/'+phase+'/'+file for phase,files in
             [('training',['predictions.json','prediction_table.json','training_geometry.json','receipt.json']),
              ('validation',['validation.json','validation_table.json','joint_scenarios.json','receipt.json'])] for file in files]
    hashes_before={p:sha(ROOT/p) for p in rels}
    protocol=read(ROOT/PROTOCOL);saved=read(BASE/'training/predictions.json');validation=read(BASE/'validation/validation.json')
    scenarios=read(BASE/'validation/joint_scenarios.json');training_geo=read(BASE/'training/training_geometry.json')
    gate=read(BASE/'prevalidation_gate.json');tr=read(BASE/'training/receipt.json');vr=read(BASE/'validation/receipt.json')
    assert protocol['status']=='frozen'
    assert gate['status']=='verified_before_first_heldout_RNA_read' and gate['RNA3_values_opened_at_this_gate'] is False
    assert gate['protocol_commit']=='faa9bb248ef3969495f54504987f63d46dcad133'
    assert gate['prediction_commit']=='d186707a74a40e14a3ed9175337796d9de25c3dd'==gate['remote_main']
    assert subprocess.run(['git','merge-base','--is-ancestor',gate['protocol_commit'],gate['prediction_commit']],cwd=ROOT).returncode==0
    committed=[]
    for commit,paths in [(gate['protocol_commit'],[PROTOCOL,RUNNER]),(gate['prediction_commit'],list(gate['sha256']))]:
        for p in paths:
            blob=subprocess.check_output(['git','show',commit+':'+p],cwd=ROOT)
            assert blob==(ROOT/p).read_bytes()
            digest=hashlib.sha256(blob).hexdigest()
            assert digest==gate['sha256'][p]
            committed.append({'commit':commit,'path':p,'sha256':digest})
    assert sha(ROOT/RUNNER)=='f57f084c85fbfae3b27888e55284ff41e10fb178d719351fe94c40d17c9c2b83'
    identity=saved['identity']
    assert validation['identity']==identity and tr['inputs']==identity
    assert vr['inputs']=={**identity,'prediction_sha256':sha(BASE/'training/predictions.json')}
    assert vr['repository']['commit']==gate['prediction_commit'] and tr['repository']['commit']==gate['protocol_commit']
    for phase,receipt in [('training',tr),('validation',vr)]:
        for name,digest in receipt['outputs'].items():assert sha(BASE/phase/name)==digest
    assert saved['exposure']['D_cell_values_accessed'] is False and validation['count_columns_accessed']==['D','G']
    keys=protocol['keys'];reference=keys.index(protocol['reference_key'])
    assert keys==saved['keys']==validation['keys'] and saved['reference_index']==reference
    assert len(keys)==len(set(keys))==32 and reference==3
    positions=[0,3,6];levels=[sorted({k[i] for k in keys}) for i in positions]
    vertices=[tuple(level.index(k[p]) for level,p in zip(levels,positions)) for k in keys]
    graph=[[j for j,y in enumerate(keys) if sum(a!=b for a,b in zip(x,y))==1] for x in keys]
    assert sum(map(len,graph))//2==112
    adjacency=np.zeros((32,32))
    for i,n in enumerate(graph):adjacency[i,n]=1
    eigenvectors_values,Q=np.linalg.eigh(7*np.eye(32)-adjacency)
    heat=(Q*np.exp(-math.log(2)/4*eigenvectors_values))@Q.T
    counts={**saved['training_counts'],**validation['validation_counts']}
    assert counts['G']==saved['training_counts']['G']
    assert set(counts)==set('BCDEFG')
    assert all(len(v)==32 and all(isinstance(n,int) and n>=0 for n in v) for v in counts.values())
    assert all(min(counts[c])>0 for c in 'EFG')
    assert all(sum(counts[c][i] for c in 'EFG')>50 for i in range(32))
    assert counts['B']!=counts['C'] and counts['D'] not in [counts['B'],counts['C']]
    auditcounts={'models':0,'vertex_predictions':0,'loss_entries':0,'geometry_fields':0,'endpoint_pairs':0,'threshold_records':0,'table_rows':0,'rectangle_identities':0}
    maxdiff={'scores':0.,'predictions':0.,'metrics':0.,'energies':0.,'scenario_field':0.,'rectangle_identity':0.}
    summaries={};dependence={};figure_arrays={}
    def audit_geometry(field,stored):
        capacity,barriers,thresholds,energies=independent_geometry(field,vertices,graph)
        assert stored['pair_count']==len(barriers)==496
        assert stored['critical_thresholds']==thresholds
        assert stored['maximum_barrier']==max(barriers.values())
        assert len(stored['pairs'])==496
        for row in stored['pairs']:
            i,j=row['i'],row['j'];assert row['capacity']==capacity[i,j] and row['barrier']==barriers[i,j]
        for o in range(4):maxdiff['energies']=max(maxdiff['energies'],check_close(energies[str(o)],stored['anova_mean_square_by_order'][str(o)]))
        auditcounts['geometry_fields']+=1;auditcounts['endpoint_pairs']+=496;auditcounts['threshold_records']+=len(thresholds)
    with localcontext() as context:
        context.prec=80
        for c in CS:
            train=score(counts,'BC','EF',c);test=score(counts,'D','G',c)
            block=saved['predictions'][c];observed=validation['results'][c]
            maxdiff['scores']=max(maxdiff['scores'],check_close(train,block['training_unanchored']),check_close(anchor(train,reference),block['training_anchored']),check_close(test,observed['test_unanchored']),check_close(anchor(test,reference),observed['test_anchored']))
            models=fit_independent(train,keys,vertices,reference,heat)
            losses={}
            for name in MODELS:
                maxdiff['predictions']=max(maxdiff['predictions'],check_close(models[name],block['models'][name]))
                metric=independent_metrics(block['models'][name],observed['test_anchored'],reference)
                for k,value in metric.items():maxdiff['metrics']=max(maxdiff['metrics'],check_close(value,observed['models'][name][k]))
                losses[name]=metric
                auditcounts['models']+=1;auditcounts['vertex_predictions']+=32;auditcounts['loss_entries']+=5
                audit_geometry(block['models'][name],training_geo[c][name])
            audit_geometry(observed['test_anchored'],observed['observed_geometry'])
            baseline_names=[n for n in MODELS if n!='pairwise_half'];best=min(losses[n]['rmse32'] for n in baseline_names)
            target=losses['pairwise_half']['rmse32'];decision=validation['summary']['all_pseudocount_decisions'][c]
            assert decision['passes_10_percent_hurdle']==(best>0 and target<=.9*best)
            assert decision['strict_improvement']==(best>0 and target<best)
            assert decision['best_baseline_names']==[n for n in baseline_names if losses[n]['rmse32']==best]
            check_close(decision['best_baseline_rmse32'],best);check_close(decision['target_rmse32'],target)
            check_close(decision['relative_improvement'],1-target/best)
            summaries[c]={'rmse32':{n:losses[n]['rmse32'] for n in MODELS},'centered_rmse':{n:losses[n]['centered_rmse'] for n in MODELS},
                          'pairwise_relative_increase':target/best-1,'passes_hurdle':decision['passes_10_percent_hurdle'],
                          'heldout_maximum_barrier':observed['observed_geometry']['maximum_barrier'],
                          'heldout_anova_mean_square_by_order':observed['observed_geometry']['anova_mean_square_by_order']}
            for r in 'BCD':
                for d in 'EFG':
                    s=scenarios[c][r+'/'+d];independent=anchor(score(counts,r,d,c),reference)
                    maxdiff['scenario_field']=max(maxdiff['scenario_field'],check_close(independent,s['field']))
                    audit_geometry(s['field'],s['geometry'])
            for r in 'CD':
                for d in 'FG':
                    residual=np.array(scenarios[c][r+'/'+d]['field'])+np.array(scenarios[c]['B/E']['field'])-np.array(scenarios[c][r+'/E']['field'])-np.array(scenarios[c]['B/'+d]['field'])
                    maxdiff['rectangle_identity']=max(maxdiff['rectangle_identity'],check_close(residual,np.zeros(32)))
                    auditcounts['rectangle_identities']+=32
            minima=min(s['geometry']['maximum_barrier'] for s in scenarios[c].values());maxima=max(s['geometry']['maximum_barrier'] for s in scenarios[c].values())
            assert validation['summary']['scenario_maximum_barrier_ranges'][c]=={'minimum':minima,'maximum':maxima}
            dependence[c]={'scenario_count':9,'four_fundamental_rectangle_identities_per_vertex':True,'span_dimension_upper_bound':5,
                           'maximum_barrier_minimum':minima,'maximum_barrier_maximum':maxima,
                           'maximum_barrier_matrix':[[scenarios[c][r+'/'+d]['geometry']['maximum_barrier'] for d in 'EFG'] for r in 'BCD']}
        for phase,table_file in [('training','prediction_table.json'),('validation','validation_table.json')]:
            rows=read(BASE/phase/table_file);assert len(rows)==672
            for row in rows:
                c=str(row['pseudocount']);i=keys.index(row['key']);name=row['model'];pred=saved['predictions'][c]['models'][name][i]
                assert row['prediction']==pred
                if phase=='validation':
                    obs=validation['results'][c]['test_anchored'][i];assert row['outcome']==obs and row['error']==pred-obs
                auditcounts['table_rows']+=1
    assert validation['summary']['primary_decision']==validation['summary']['all_pseudocount_decisions']['0.5']
    assert all(not summaries[c]['passes_hurdle'] for c in CS)
    assert validation['summary']['primary_success_and_sensitivity_direction'] is False
    primary=validation['results']['0.5'];figure_arrays={'model_order':MODELS,'primary_rmse32':[primary['models'][n]['rmse32'] for n in MODELS],
        'ten_percent_hurdle':.9*min(primary['models'][n]['rmse32'] for n in MODELS if n!='pairwise_half'),
        'training_anova_mean_square':training_geo['0.5']['training_mean']['anova_mean_square_by_order'],
        'validation_anova_mean_square':primary['observed_geometry']['anova_mean_square_by_order'],
        'scenario_row_order':['B','C','D'],'scenario_column_order':['E','F','G'],
        'scenario_maximum_barrier_matrix':dependence['0.5']['maximum_barrier_matrix']}
    hashes_after={p:sha(ROOT/p) for p in rels};assert hashes_before==hashes_after
    result={'status':'passed','scope':'Independent audit of frozen32-count/prediction vectors only, no workbook opened',
            'pseudocounts':CS,'primary_result':'Predeclared pairwise predictor failed; all sensitivity settings also failed',
            'summaries':summaries,'scenario_dependence':dependence,'figure_arrays':figure_arrays,
            'checked':auditcounts,'largest_absolute_difference':maxdiff,'absolute_numeric_tolerance':TOL,
            'decimal_logarithm_precision':80,'committed_blob_checks':committed,
            'chronology_limit':'Committed content and recorded prevalidation gate are consistent; this audit did not independently observe the historical first RNA3 read.',
            'input_hashes_before':hashes_before,'input_hashes_after':hashes_after,
            'script_sha256':sha(Path(__file__)),'observed_workbook_opened':False,
            'started_at_utc':started_at,'runtime_seconds':time.perf_counter()-started,'python':sys.version,'numpy':np.__version__}
    out.mkdir(parents=True)
    (out/'audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    code=Path(__file__).relative_to(ROOT).as_posix();output=(out/'audit.json').relative_to(ROOT).as_posix()
    manifest={'schema_version':2,'claim_id':'GILL_FROZEN_PILOT_INDEPENDENT_RESULT_AUDIT',
        'repository':{'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                      'dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True))},
        'command':[sys.executable,*sys.argv],
        'environment':{'software':[{'name':'CPython','version':sys.version.split()[0]},{'name':'NumPy','version':np.__version__}],
                       'hardware':sys.platform},
        'mathematics':{'assertion_tested':'All frozen saved32-vector formulas, models, losses, graph geometries, table values and recorded freeze identities reproduce independently',
                       'coefficient_domain':'Integer counts;80-digit Decimal logarithms;float64 saved fields;exact finite threshold connectivity',
                       'conventions':'Fixed32genotypes;no observed workbook access;all3predeclared pseudocounts;all7frozen models;no retuning',
                       'inputs':[code,*rels],'bounds':auditcounts,
                       'non_claims':['No calibrated biological inference or new-genotype generalization','The historical first RNA3 read was not independently observed','No claim of signal/error separation']},
        'randomness':{'used':False},
        'run':{'started_at':started_at,'runtime_seconds':result['runtime_seconds'],'status':'completed','exit_status':0,
               'resource_limits':{'memory_bytes':None,'cpu_seconds':None,'max_cores':None,'max_threads':None}},
        'input_artifacts':[{'path':p,'sha256':hashes_before[p],'sha256_after':hashes_after[p]} for p in rels]
                          +[{'path':code,'sha256':result['script_sha256'],'sha256_after':sha(Path(__file__))}],
        'outputs':[{'path':output,'sha256':sha(out/'audit.json'),'kind':'result'}],
        'checks':[{'name':'source-and-saved-inputs-unchanged','passed':True},{'name':'all-predeclared-result-checks','passed':True}],
        'result':'completed','declared_results':[output],
        'residual_risks':['Floating-point comparisons use absolute tolerance2e-12, with measured discrepancies recorded','No hard resource limits; finite deterministic audit']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'checked':auditcounts,'differences':maxdiff,'primary':summaries['0.5'],'runtime_seconds':result['runtime_seconds']},indent=2))

if __name__=='__main__':main()
