"""Frozen hypothetical missingness calculations. Metadata inputs only."""
if not __debug__:
    raise SystemExit('Refuse disabled assertions before input access')
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import numpy as np

DESIGN_SHA='90a3e1cdeb425cbaf482df4de0ecac3c701510b0fa95e6836918428752ba5363'
ARMS=('development','reserved_validation')
STAGES=('descriptive','known_corner_prediction','prediction_evaluation')

def check(condition,message):
    if not condition:raise RuntimeError(message)

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def save(path,data):
    with path.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(data,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')

def prepare(split,mask):
    units=sorted(split['units'],key=lambda u:u['square_index'])
    check(len(units)==47,'All47 squares required')
    groups=sorted([','.join(sorted(g,key=int)) for g in split['seed_groups']],key=lambda x:tuple(map(int,x.split(','))))
    by_seed={seed:','.join(sorted(g,key=int)) for g in split['seed_groups'] for seed in g}
    ids=sorted({i for u in units for i in u['coordinate_order_00_01_10_11']})
    check(len(ids)==181 and len(groups)==32,'Metadata dimensions changed')
    idx={v:i for i,v in enumerate(ids)};global_mask={v for a in ARMS for v in mask['arms'][a]['masked_IDs']}
    check(len(global_mask)==46,'Global mask changed')
    mask_by_index={u['square_index']:u for u in mask['units']}
    id_group={};structural=[];known=[]
    for u in units:
        corners=u['coordinate_order_00_01_10_11'];g=by_seed[u['seed']]
        for v in corners:
            check(id_group.get(v,g)==g,'Shared ID crosses groups');id_group[v]=g
        eligible=len(set(corners)&global_mask)==1
        check(mask_by_index[u['square_index']]['prediction_eligible_from_metadata']==eligible,'Mask eligibility mismatch')
        structural.append(eligible)
        # A dummy 3-entry index for an excluded case is masked out after evaluation.
        known.append([idx[v] for v in corners if v not in global_mask] if eligible else [0,0,0])
    return {'units':units,'ids':ids,'groups':groups,'global_mask':global_mask,'by_seed':by_seed,
        'corners':np.array([[idx[v] for v in u['coordinate_order_00_01_10_11']] for u in units]),
        'known':np.array(known),'structural':np.array(structural,dtype=bool),
        'id_group_index':np.array([groups.index(id_group[v]) for v in ids]),
        'masked_index':np.array([idx[v] for v in sorted(global_mask)]),
        'arm_positions':{a:[i for i,u in enumerate(units) if u['arm']==a] for a in ARMS}}

def availability(meta,family,p,rng,n):
    if family=='A':return np.all(rng.random((n,len(meta['ids']),7))>=p,axis=2)
    if family=='B':return (rng.random((n,len(meta['groups'])))>=p)[:,meta['id_group_index']]
    if family=='C':
        present=np.ones((n,len(meta['ids'])),dtype=bool)
        present[:,meta['masked_index']]=rng.random((n,len(meta['masked_index'])))>=p
        return present
    raise ValueError('Unknown frozen family')

def vector_metrics(meta,present):
    desc=np.all(present[:,meta['corners']],axis=2)
    pred=np.all(present[:,meta['known']],axis=2)&meta['structural'][None,:]
    evaluated=desc&meta['structural'][None,:];unscorable=pred&~evaluated
    check(not np.any(evaluated&~pred),'Evaluation without a prediction')
    result={}
    for arm in ARMS:
        positions=meta['arm_positions'][arm];units=[meta['units'][i] for i in positions]
        groups=sorted({meta['by_seed'][u['seed']] for u in units})
        targets=sorted({u['withheld_ID'] for u in units})
        arm_result={}
        for stage,eligible in zip(STAGES,(desc,pred,evaluated)):
            local=eligible[:,positions]
            count=local.sum(axis=1)
            group_count=np.sum([np.any(local[:,[j for j,u in enumerate(units) if meta['by_seed'][u['seed']]==g]],axis=1) for g in groups],axis=0)
            target_count=np.sum([np.any(local[:,[j for j,u in enumerate(units) if u['withheld_ID']==t]],axis=1) for t in targets],axis=0)
            check(np.all(target_count<=count),'Target count exceeds case count')
            arm_result[stage+'_squares']=count;arm_result[stage+'_groups']=group_count
            arm_result[stage+'_distinct_targets']=target_count
            arm_result[stage+'_coverage_pass']=(count>=10)&(group_count>=8)
        local=unscorable[:,positions]
        arm_result['predicted_but_unscorable_cases']=local.sum(axis=1)
        arm_result['distinct_unscorable_targets']=np.sum([np.any(local[:,[j for j,u in enumerate(units) if u['withheld_ID']==t]],axis=1) for t in targets],axis=0)
        check(np.all(arm_result['known_corner_prediction_squares']==arm_result['prediction_evaluation_squares']+arm_result['predicted_but_unscorable_cases']),'Prediction accounting fails')
        result[arm]=arm_result
    return result

def scalar_reference(meta,present):
    """Distinct Python set/loop implementation, no numerical indexing formulas."""
    available={identifier for identifier,is_present in zip(meta['ids'],present) if bool(is_present)}
    result={}
    for arm in ARMS:
        result[arm]={};lists={k:[] for k in STAGES};missing=[]
        for u in meta['units']:
            if u['arm']!=arm:continue
            corners=set(u['coordinate_order_00_01_10_11']);hits=corners&meta['global_mask']
            is_descriptive=corners<=available
            is_prediction=len(hits)==1 and corners-hits<=available
            is_evaluation=len(hits)==1 and is_descriptive
            if is_descriptive:lists['descriptive'].append(u)
            if is_prediction:lists['known_corner_prediction'].append(u)
            if is_evaluation:lists['prediction_evaluation'].append(u)
            if is_prediction and not is_evaluation:missing.append(u)
        for stage,rows in lists.items():
            n=len(rows);g=len({meta['by_seed'][u['seed']] for u in rows})
            result[arm][stage+'_squares']=n;result[arm][stage+'_groups']=g
            result[arm][stage+'_distinct_targets']=len({u['withheld_ID'] for u in rows})
            result[arm][stage+'_coverage_pass']=n>=10 and g>=8
        result[arm]['predicted_but_unscorable_cases']=len(missing)
        result[arm]['distinct_unscorable_targets']=len({u['withheld_ID'] for u in missing})
    return result

def verify_reference(meta,present,metrics,limit):
    comparisons=0
    for i in range(min(limit,len(present))):
        independent=scalar_reference(meta,present[i])
        for arm in ARMS:
            for key,value in independent[arm].items():
                check(value==metrics[arm][key][i],f'Scalar disagreement {arm}/{key}/{i}')
                comparisons+=1
    return comparisons

def controls(meta):
    out=[]
    for family in ('A','B','C'):
        for p in (0,1):
            # Deterministic boundary inputs, no stochastic draws required.
            present=np.ones((1,len(meta['ids'])),dtype=bool)
            if p==1:
                if family=='C':present[:,meta['masked_index']]=False
                else:present[:]=False
            values=vector_metrics(meta,present)
            verify_reference(meta,present,values,1)
            for arm in ARMS:
                rows=[u for u in meta['units'] if u['arm']==arm]
                n=len(rows);eligible=sum(len(set(u['coordinate_order_00_01_10_11'])&meta['global_mask'])==1 for u in rows)
                check(values[arm]['descriptive_squares'][0]==(n if p==0 else 0),'Boundary descriptive failure')
                check(values[arm]['prediction_evaluation_squares'][0]==(eligible if p==0 else 0),'Boundary evaluation failure')
                check(values[arm]['known_corner_prediction_squares'][0]==(eligible if p==0 or family=='C' else 0),'Boundary prediction failure')
                check(values[arm]['predicted_but_unscorable_cases'][0]==(eligible if p==1 and family=='C' else 0),'Boundary unscorable failure')
            out.append({'family':family,'p':p,'metrics':{a:{k:v[0].item() for k,v in d.items()} for a,d in values.items()}})
    return out

def rate(events,condition=None):
    if condition is not None:events=events[condition]
    n=len(events);count=int(np.sum(events));estimate=count/n if n else None
    return {'successes':count,'denominator':n,'rate':estimate,
        'Monte_Carlo_standard_error':math.sqrt(estimate*(1-estimate)/n) if n else None}

def summarize(values,probs):
    result={}
    for key,value in values.items():
        if key.endswith('_coverage_pass'):result[key]=rate(value)
        else:
            sd=float(np.std(value,ddof=1))
            result[key]={'mean':float(np.mean(value)),'sample_standard_deviation':sd,
                'Monte_Carlo_standard_error_of_mean':sd/math.sqrt(len(value)),
                'quantiles':{str(q):int(v) for q,v in zip(probs,np.quantile(value,probs,method='inverted_cdf'))}}
    result['evaluation_given_descriptive_coverage']=rate(values['prediction_evaluation_coverage_pass'],values['descriptive_coverage_pass'])
    result['evaluation_given_known_corner_coverage']=rate(values['prediction_evaluation_coverage_pass'],values['known_corner_prediction_coverage_pass'])
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--design',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();check(not args.output.exists(),'Output directory exists')
    check(sha(args.design)==DESIGN_SHA,'Frozen design changed')
    design=json.loads(args.design.read_text())
    check(platform.python_version()==design['runtime']['python'] and np.__version__==design['runtime']['numpy'],'Frozen runtime mismatch')
    inputs={}
    for role,entry in design['inputs'].items():
        path=(args.design.parent/entry['path']).resolve();check(sha(path)==entry['sha256'],'Input changed: '+role)
        inputs[role]=json.loads(path.read_text())
    meta=prepare(inputs['split'],inputs['mask']);exact=controls(meta)
    args.output.mkdir(parents=True,exist_ok=False)
    freeze={'status':'frozen_before_random_draws','design_sha256':DESIGN_SHA,'source_sha256':sha(Path(__file__)),
        'input_hashes':{r:e['sha256'] for r,e in design['inputs'].items()},'python':sys.version,'numpy':np.__version__,
        'ordered_genotype_IDs':meta['ids'],'ordered_group_IDs':meta['groups'],'sorted_masked_IDs':sorted(meta['global_mask'])}
    save(args.output/'execution_freeze.json',freeze)
    save(args.output/'exact_controls.json',exact)
    results=[];references=0;flat=[]
    for index,scenario in enumerate(design['scenarios']):
        rng=np.random.Generator(np.random.PCG64(scenario['seed']))
        blocks={a:{} for a in ARMS};remaining=design['draws_per_scenario'];first=True
        while remaining:
            n=min(design['runtime']['batch_size'],remaining)
            present=availability(meta,scenario['family'],scenario['p'],rng,n);current=vector_metrics(meta,present)
            if first:references+=verify_reference(meta,present,current,64);first=False
            for arm in ARMS:
                for key,v in current[arm].items():blocks[arm].setdefault(key,[]).append(v)
            remaining-=n
        values={a:{k:np.concatenate(v) for k,v in d.items()} for a,d in blocks.items()}
        filename=f'scenario_{index+1:02d}_counts.npz'
        np.savez_compressed(args.output/filename,**{a+'__'+k:v for a,d in values.items() for k,v in d.items()})
        summary={a:summarize(d,design['reporting']['quantile_probabilities']) for a,d in values.items()}
        joint={stage:rate(values[ARMS[0]][stage+'_coverage_pass']&values[ARMS[1]][stage+'_coverage_pass']) for stage in STAGES}
        results.append({**scenario,'draws':design['draws_per_scenario'],'arms':summary,'both_arms_pass':joint,
            'count_array_file':filename,'count_array_sha256':sha(args.output/filename)})
        for arm,d in summary.items():
            flat.append({'family':scenario['family'],'p':scenario['p'],'arm':arm,
                'descriptive_pass_rate':d['descriptive_coverage_pass']['rate'],
                'prediction_stage_pass_rate':d['known_corner_prediction_coverage_pass']['rate'],
                'prediction_evaluation_pass_rate':d['prediction_evaluation_coverage_pass']['rate'],
                'evaluation_given_descriptive':d['evaluation_given_descriptive_coverage']['rate'],
                'evaluation_given_prediction_stage':d['evaluation_given_known_corner_coverage']['rate'],
                'mean_descriptive_squares':d['descriptive_squares']['mean'],
                'mean_evaluated_squares':d['prediction_evaluation_squares']['mean'],
                'mean_evaluated_distinct_targets':d['prediction_evaluation_distinct_targets']['mean'],
                'mean_predicted_unscorable_cases':d['predicted_but_unscorable_cases']['mean']})
    result={'status':'all_frozen_scenarios_completed','design_sha256':DESIGN_SHA,'source_sha256':sha(Path(__file__)),
        'scenarios':results,'total_synthetic_draws':sum(s['draws'] for s in results),
        'reference_draws':64*len(results),'reference_scalar_comparisons':references,'deterministic_boundary_controls':6,
        'all_stage_containment_checks_passed':True,'scope':'Hypothetical coverage only; no biological uncertainty or power estimate.',
        'real_phenotype_values_read':0}
    save(args.output/'results.json',result)
    with (args.output/'summary.csv').open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    save(args.output/'receipt.json',{'status':result['status'],'files':{p.name:sha(p) for p in args.output.iterdir() if p.is_file()},
        'design_sha256':DESIGN_SHA,'source_sha256':sha(Path(__file__)),'real_phenotype_values_read':0})
    print(json.dumps({'status':result['status'],'synthetic_draws':result['total_synthetic_draws'],
        'reference_comparisons':references,'results_sha256':sha(args.output/'results.json')}))

if __name__=='__main__':main()
