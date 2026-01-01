"""Hash/arm guarded finite replay from explicitly authorized selected exports.

No archive reader, network, Git query, source discovery or phenotype path defaults.
"""
if not __debug__:
    raise SystemExit('Disabled assertions refused before imports or input/output')
import argparse
import csv
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True

REFERENCE_SHA='64a5e9040a77d81dc779583d1d3ea5894a5aebd5718dbf3466281803f48ce375'
PROTEIN=['clean.lin.prot.repOne','clean.lin.prot.repTwo','clean.lin.prot.repThree','clean.lin.prot.repFour']
GROWTH=['growth.R60one','growth.R60two','growth.R60three']
REPS=PROTEIN+GROWTH
MEANS=['clean.lin.prot.mean','growth.mean.R60']
METHODS={'additive':'additive','neighbor_mean':'adjacent_mean','three_mean':'three_corner_mean'}
LABELS={'joint_connected':'jointly_connected','joint_only':'joint_only_obstruction','individual_disconnected':'single_trait_disconnection','ineligible_endpoint':'ineligible_endpoints'}
ROLES={'evaluation','unsealed_inputs','predictions','prediction_seal','known_inputs','split','components','edges','adapter','reference'}
BASE_ROLES=['protocol','split','mask','components','edges','archive','graph_rules','masking','runner']
TOLERANCE=2e-12

def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def unique_object(pairs):
    out={}
    for k,v in pairs:
        if k in out:raise ValueError('Duplicate JSON key: '+k)
        out[k]=v
    return out
def load(path):
    with Path(path).open(encoding='utf-8') as stream:
        return json.load(stream,object_pairs_hook=unique_object,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON')))
def write(path,x):
    with Path(path).open('x',encoding='utf-8') as stream:json.dump(x,stream,indent=2,allow_nan=False);stream.write('\n')

class Check:
    def __init__(self):self.exact=0;self.numeric=0
    def equal(self,actual,expected,where):
        if where.rsplit('.',1)[-1] in ('thresholds','point_tolerance','interval_tolerance','tolerance'):
            if actual!=expected:raise ValueError('Exact threshold/tolerance mismatch: '+where)
            self.exact+=1;return
        if type(expected) is dict:
            if not isinstance(actual,dict) or set(actual)!=set(expected):raise ValueError('Dictionary keys differ: '+where)
            for k,v in expected.items():self.equal(actual[k],v,where+'.'+str(k))
        elif type(expected) is list:
            if not isinstance(actual,list) or len(actual)!=len(expected):raise ValueError('List shape differs: '+where)
            for k,v in enumerate(expected):self.equal(actual[k],v,where+'['+str(k)+']')
        elif isinstance(expected,float):
            if isinstance(actual,bool) or not isinstance(actual,(int,float)) or not math.isfinite(actual) or not math.isclose(actual,expected,rel_tol=TOLERANCE,abs_tol=TOLERANCE):
                raise ValueError('Numerical mismatch: '+where+' actual='+repr(actual)+' expected='+repr(expected))
            self.numeric+=1
        else:
            if actual!=expected or isinstance(expected,bool)!=isinstance(actual,bool):raise ValueError('Exact mismatch: '+where+' actual='+repr(actual)+' expected='+repr(expected))
            self.exact+=1
    def fields(self,actual,expected,where):
        for k,v in expected.items():
            if k not in actual:raise ValueError('Missing field: '+where+'.'+k)
            self.equal(actual[k],v,where+'.'+k)

def authorize(path,root):
    auth=load(path)
    if auth.get('schema')!='Cambray_independent_result_v1' or auth.get('execution_authorized') is not True:
        raise ValueError('Explicit comparison authorization required')
    if auth.get('data_scope') not in ('synthetic','explicitly_authorized_observed'):raise ValueError('Undeclared data scope')
    if auth.get('arm') not in ('development','reserved_validation'):raise ValueError('Invalid arm')
    if set(auth.get('inputs',{}))!=ROLES:raise ValueError('Require exactly the declared input roles')
    paths={};hashes={}
    for role,entry in auth['inputs'].items():
        p=Path(entry['path']);p=(p if p.is_absolute() else root/p).resolve()
        if sha(p)!=entry['sha256']:raise ValueError('Hash mismatch: '+role)
        paths[role]=p;hashes[role]=entry['sha256']
    if paths['adapter']!=Path(__file__).resolve():raise ValueError('Pinned adapter is not executing adapter')
    if hashes['reference']!=REFERENCE_SHA:raise ValueError('Independent reference changed')
    spec=importlib.util.spec_from_file_location('frozen_independent_reference',paths['reference'])
    ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
    return auth,paths,hashes,ref

def validate_row(row):
    if set(row)!=set(REPS+MEANS):raise ValueError('Unexpected selected trait columns')
    for cell in row.values():
        if set(cell)!={'value','missing_reason','raw'} or not isinstance(cell['raw'],str):raise ValueError('Malformed selected cell')
        value=cell['value']
        if value is not None and (isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value)):
            raise ValueError('Nonfinite selected value')
        if value is None:
            if cell['missing_reason'] is None:raise ValueError('Missing value lacks reason')
        elif cell['missing_reason'] is not None or float(cell['raw'])!=value:raise ValueError('Raw numeric/value discrepancy')

def mean(row,columns):
    values=[row[c]['value'] for c in columns]
    return None if any(v is None for v in values) else math.fsum(values)/len(values)
def point(row):return [mean(row,c) for c in (PROTEIN,GROWTH)]
def ranges(row):
    return [None if any(row[c]['value'] is None for c in columns) else [min(row[c]['value'] for c in columns),max(row[c]['value'] for c in columns)] for columns in (PROTEIN,GROWTH)]
def missing(rows,ids):return [{'ID':i,'column':c,'reason':rows[i][c]['missing_reason']} for i in ids for c in REPS if rows[i][c]['value'] is None]
def tolerances(field):return [1e-10*max([1]+[abs(row[k]) for row in field.values() if row[k] is not None]) for k in (0,1)]
def range_tolerance(field):return [1e-10*max([1]+[abs(x) for row in field.values() if row[k] is not None for x in row[k]]) for k in (0,1)]

def classify(ref,field,edges,endpoints,tol=(0,0)):
    r=ref.classify(field,edges,*endpoints,tolerances=tol)
    return {'classification':LABELS[r['status']],'thresholds':r['thresholds'],'individual_connectivity':r['individual'],'joint_connectivity':r['joint']}
def bounds(ref,field,edges,endpoints,tol=(0,0)):
    tau=[min(field[v][k][0] for v in endpoints) for k in (0,1)]
    # Shift thresholds, not scores, to preserve the protocol's exact f>=tau-delta operation.
    r=ref.bounds(field,edges,*endpoints,[tau[k]-tol[k] for k in (0,1)])
    if r['status']=='ineligible_endpoint':return {'status':'ineligible_definite_endpoints','thresholds':tau}
    if r['joint_pessimistic']:status='certified_jointly_connected'
    elif not all(r['individual_optimistic']):status='certified_single_trait_disconnection'
    elif r['joint_only_certificate']:status='certified_joint_only_obstruction'
    else:status='inconclusive'
    return dict(status=status,thresholds=tau,**{k:r[k] for k in ('individual_pessimistic','individual_optimistic','joint_pessimistic','joint_optimistic')})
def field_check(ref,field,edges,endpoints,tol=(0,0)):
    complete=all(x is not None for row in field.values() for x in row)
    intervals={v:[None if x is None else [x,x] for x in row] for v,row in field.items()}
    return {'complete':complete,'point':classify(ref,field,edges,endpoints,tol) if complete else None,'bounds':bounds(ref,intervals,edges,endpoints,tol)}
def rational(value):
    if isinstance(value,dict) and set(value)=={'numerator','denominator'}:return F(value['numerator'],value['denominator'])
    raise ValueError('Malformed exact geometry coordinate')

def check_geometry(check,ref,observed,field,endpoints):
    exact={v:[F(x) for x in row] for v,row in field.items()};caps=[min(exact[v][k] for v in endpoints) for k in (0,1)]
    axes=[]
    for k in (0,1):
        values=sorted({row[k] for row in exact.values() if row[k]<=caps[k]})
        axis=[('open',None,values[0],values[0]-1)]
        for n,value in enumerate(values):
            axis.append(('atom',value,value,value))
            if n+1<len(values):axis.append(('open',value,values[n+1],(value+values[n+1])/2))
        axes.append(axis)
    if len(observed['axes'])!=2:raise ValueError('Geometry axis count')
    for k,axis in enumerate(axes):
        check.equal(len(observed['axes'][k]),len(axis),'geometry.axis_length')
        for o,(kind,lo,hi,rep) in zip(observed['axes'][k],axis):
            check.equal(o['kind'],kind,'geometry.kind')
            if (None if o['lower'] is None else rational(o['lower']))!=lo or rational(o['upper'])!=hi or rational(o['representative'])!=rep:
                raise ValueError('Exact geometry boundary mismatch')
    if list(map(rational,observed['endpoint_caps']))!=caps:raise ValueError('Geometry endpoint caps')
    expected=[]
    for i,a in enumerate(axes[0]):
        for j,b in enumerate(axes[1]):
            # The frozen reference's point type guard accepts int/float only.
            # Keep exact Fraction comparisons here, then use its independent closure.
            allowed=[{v for v,row in exact.items() if row[k]>=(a[3],b[3])[k]} for k in (0,1)]
            individual=[ref.closure_connected(list(exact),ref.SQUARE_EDGES,subset,*endpoints) for subset in allowed]
            joint=ref.closure_connected(list(exact),ref.SQUARE_EDGES,allowed[0]&allowed[1],*endpoints)
            label='joint_connected' if joint else 'joint_only' if all(individual) else 'individual_disconnected'
            expected.append({'axis_indices':[i,j],'classification':LABELS[label]})
    check.equal(observed['cells'],expected,'geometry.cells')
    return len(expected)

def stats():return dict(eligible_squares=0,eligible_diagonals=0,joint_only_obstructions=0,preparation_consistent_diagonals=0,range_certified_diagonals=0)
def fraction(s):s['joint_only_fraction']=s['joint_only_obstructions']/s['eligible_diagonals'] if s['eligible_diagonals'] else None

def compare(paths,hashes,arm,ref,check):
    # All provided hashes are checked by authorize before any selected outcome is decoded.
    seal=load(paths['prediction_seal'])
    check.fields(seal,{'schema':'Cambray_prediction_seal_v1','arm':arm,'status':'predictions_sealed_before_target_access',
        'target_finiteness_examined':False,'author_means_examined':False,'component_outcomes_examined':False},'seal')
    check.equal(set(seal['files']),{'known_inputs.json','predictions.json'},'sealed filenames')
    check.equal(seal['files']['predictions.json'],hashes['predictions'],'sealed predictions hash')
    check.equal(seal['files']['known_inputs.json'],hashes['known_inputs'],'sealed known inputs hash')
    for role in ('split','components','edges'):check.equal(seal['input_hashes'][role],hashes[role],'seal metadata '+role)
    evaluation=load(paths['evaluation']);unsealed=load(paths['unsealed_inputs']);predictions=load(paths['predictions']);known=load(paths['known_inputs']);split=load(paths['split'])
    for name,obj in (('evaluation',evaluation),('unsealed',unsealed),('predictions',predictions),('known',known)):check.equal(obj['arm'],arm,name+' arm')
    check.equal(evaluation['prediction_seal_sha256'],hashes['prediction_seal'],'evaluation seal')
    check.equal(predictions['target_finiteness_examined'],False,'prediction no finiteness')
    check.equal(predictions['all_source_and_metadata_hashes'],seal['input_hashes'],'prediction inputs')
    for role in BASE_ROLES:check.equal(evaluation['source_hashes'][role],seal['input_hashes'][role],'evaluation inputs '+role)
    components=[json.loads(line,object_pairs_hook=unique_object) for line in paths['components'].read_text(encoding='utf-8').splitlines()]
    units=[u for u in split['units'] if u['arm']==arm]
    if not units or len({u['square_index'] for u in units})!=len(units):raise ValueError('Empty or duplicate unit inventory')
    group={seed:','.join(g) for g in split['seed_groups'] for seed in g}
    cids={u['component_index']:sorted(components[u['component_index']]['IDs']) for u in units}
    all_ids={i for ids in cids.values() for i in ids};rows=unsealed['rows']
    check.equal(set(rows),all_ids,'unsealed ID whitelist')
    check.fields(unsealed['reader'],{'selected_records':len(all_ids),'decoded_fields':REPS+MEANS,'other_row_phenotype_values_decoded':0},'unsealed reader')
    for row in rows.values():validate_row(row)
    mask={u['withheld_ID'] for u in units};check.equal(predictions['masked_IDs'],sorted(mask),'union mask')
    structural={u['square_index']:len(set(u['coordinate_order_00_01_10_11'])&mask)==1 for u in units}
    visible={v for u in units if structural[u['square_index']] for v in u['coordinate_order_00_01_10_11'] if v not in mask}
    check.equal(set(known['rows']),visible,'known ID whitelist')
    check.fields(known['reader'],{'selected_records':len(visible),'decoded_fields':REPS,'other_row_phenotype_values_decoded':0},'known reader')
    for v in visible:check.equal(known['rows'][v],{c:rows[v][c] for c in REPS},'known rows '+v)
    edges=[]
    with paths['edges'].open(encoding='utf-8',newline='') as stream:
        for edge in csv.DictReader(stream,delimiter='\t'):
            a,b=edge['left_ID'],edge['right_ID']
            if a in all_ids and b in all_ids:edges.append((a,b))
    coverage=[];mean_checks=[];holds=[]
    for v in sorted(rows):
        for k,cols in enumerate((PROTEIN,GROWTH)):
            value=mean(rows[v],cols);author=rows[v][MEANS[k]]['value']
            if value is None:status='not_comparable_missing_replicates'
            elif author is None:status='missing_author_mean';holds.append([v,('protein','growth')[k],status])
            elif math.isclose(value,author,rel_tol=1e-8,abs_tol=1e-10):status='compatible'
            else:status='mean_discrepancy';holds.append([v,('protein','growth')[k],status])
            mean_checks.append(dict(ID=v,trait=('protein','growth')[k],computed=value,author=author,
                difference=None if value is None or author is None else value-author,
                tolerance=None if value is None or author is None else max(1e-10,1e-8*max(abs(value),abs(author))),status=status))
    check.equal(evaluation['mean_diagnostics'],mean_checks,'mean diagnostics')
    for u in units:
        ids=u['coordinate_order_00_01_10_11']
        if len(set(ids))!=4 or u['withheld_ID'] not in ids or not set(ids)<=set(cids[u['component_index']]):raise ValueError('Invalid square membership')
        miss=missing(rows,ids)
        coverage.append(dict(square_index=u['square_index'],seed=u['seed'],group=group[u['seed']],eligible=not miss,missing_columns=miss,
            all_diagonals_retained=u['diagonals'],prediction_metadata_eligible=structural[u['square_index']]))
    check.equal(evaluation['coverage'],coverage,'coverage')
    if holds:
        check.fields(evaluation,{'status':'held_for_author_mean_provenance','interpretation_permitted':False,'holds':holds},'held evaluation')
        if any(k in evaluation for k in ('descriptions','prediction_evaluation','primary_summary')):raise ValueError('Scientific summaries present after mean hold')
        return {'status':'verified_source_hold','mean_checks':len(mean_checks),'holds':holds,'coverage':coverage}
    check.fields(evaluation,{'status':'arm_evaluation_complete','interpretation_permitted':True},'evaluation status')
    predmap={p['square_index']:p for p in predictions['records']}
    check.equal(len(predmap),len(predictions['records']),'duplicate prediction records')
    check.equal(list(predmap),[u['square_index'] for u in units],'prediction unit order')
    expected_cases=[];case_map={};scored_records=[]
    for u in units:
        p=predmap[u['square_index']];ids=u['coordinate_order_00_01_10_11'];where='prediction '+str(u['square_index'])
        check.fields(p,{'seed':u['seed'],'group':group[u['seed']],'withheld_ID':u['withheld_ID']},where)
        if not structural[u['square_index']]:
            check.fields(p,{'status':'excluded_shared_corner_mask','reason':'multiple_globally_masked_corners'},where);status=p['status']
        elif missing(rows,[v for v in ids if v not in mask]):
            check.fields(p,{'status':'missing_known_replicates','missing_columns':missing(rows,[v for v in ids if v not in mask])},where);status=p['status']
        else:
            target=ids.index(u['withheld_ID']);endpoints=(target^1,target^2);known_field={j:point(rows[v]) for j,v in enumerate(ids) if j!=target}
            estimate=ref.local_predictions(known_field,target)
            pc={m:ref.classify({j:(value if j==target else known_field[j]) for j in range(4)},ref.SQUARE_EDGES,*endpoints)['status']=='joint_only' for m,value in estimate.items()}
            check.fields(p,{'status':'predicted','visible_IDs':[v for v in ids if v not in mask],
                'known_corner_means':{v:point(rows[v]) for v in ids if v not in mask},
                'predictions':{METHODS[m]:v for m,v in estimate.items()},
                'classification_endpoint_IDs':[ids[j] for j in endpoints]},where)
            for m,primary in METHODS.items():
                result=ref.classify({j:(estimate[m] if j==target else known_field[j]) for j in range(4)},ref.SQUARE_EDGES,*endpoints)
                check.equal(p['classifications'][primary],{'classification':LABELS[result['status']],'joint_only_obstruction':pc[m]},where+' class '+m)
            status='unscorable_missing_target_or_square' if missing(rows,ids) else 'scored'
            if status=='scored':
                true=ref.classify({j:point(rows[v]) for j,v in enumerate(ids)},ref.SQUARE_EDGES,*endpoints)['status']=='joint_only'
                case=dict(group=group[u['seed']],truth=point(rows[u['withheld_ID']]),predictions=estimate,true_class=true,predicted_class=pc)
                expected_cases.append(case);case_map[u['square_index']]=case
        scored_records.append((u,p,status))
    score=evaluation['prediction_evaluation'];check.equal(len(score['records']),len(units),'score record count')
    for actual,(u,p,status) in zip(score['records'],scored_records):
        check.fields(actual,{'square_index':u['square_index'],'seed':u['seed'],'group':group[u['seed']],'status':status},'score record')
        if status in ('excluded_shared_corner_mask','missing_known_replicates'):
            check.equal(actual['reason'],p.get('reason',p.get('missing_columns')),'unscored reason')
        if status=='unscorable_missing_target_or_square':check.equal(actual['missing_columns'],missing(rows,u['coordinate_order_00_01_10_11']),'missing target record')
        if status=='scored':
            case=case_map[u['square_index']]
            observed=case['truth'];check.fields(actual,{'observed_target':observed,'true_joint_only_obstruction':case['true_class']},'score truth')
            for method,primary in METHODS.items():
                errors=[a-b for a,b in zip(p['predictions'][primary],observed)]
                check.fields(actual['methods'][primary],{'predicted_target':p['predictions'][primary],'signed_error':errors,
                    'absolute_error':list(map(abs,errors)),'squared_error':[x*x for x in errors],
                    'predicted_joint_only_obstruction':case['predicted_class'][method],
                    'classification_error':int(case['predicted_class'][method]!=case['true_class'])},'score error '+primary)
    metrics=ref.grouped_metrics(expected_cases)
    check.fields(score,{'scored_squares':metrics['n_cases'],'related_seed_groups':metrics['n_groups'],'both_true_classes':metrics['both_classes'],'discrimination_screen_passed':metrics['success']},'score summary')
    for m,primary in METHODS.items():
        expected=metrics['metrics'][m];actual=score['summaries'][primary]
        grouped={g:[c for c in expected_cases if c['group']==g] for g in sorted({c['group'] for c in expected_cases})}
        for k,trait in enumerate(('protein','growth')):
            ga={g:math.fsum(abs(c['predictions'][m][k]-c['truth'][k]) for c in cs)/len(cs) for g,cs in grouped.items()}
            gs={g:math.fsum((c['predictions'][m][k]-c['truth'][k])**2 for c in cs)/len(cs) for g,cs in grouped.items()}
            check.fields(actual['continuous'][trait],{'MAE':expected['MAE'][k],'RMSE':expected['RMSE'][k],'group_MAE':ga,'group_MSE':gs},'metric '+primary+' '+trait)
        check.fields(actual,{'pooled_confusion':expected['confusion'],'pooled_balanced_accuracy':expected['balanced_accuracy']},'metric class '+primary)
        gc={g:sum(c['predicted_class'][m]!=c['true_class'] for c in cs)/len(cs) for g,cs in grouped.items()}
        check.equal(actual['group_balanced_classification_error'],{'value':expected['misclassification'],'groups':gc},'misclassification '+primary)
    expected_comparisons=[];reasons=[]
    if not metrics['coverage_pass']:reasons.append('prediction_coverage_below_10_squares_8_groups')
    if not metrics['both_classes']:reasons.append('both_true_classes_not_observed')
    for independent,primary in list(METHODS.items())[1:]:
        for k,trait in enumerate(('protein','growth')):
            candidate=metrics['metrics']['additive']['MAE'][k];baseline=metrics['metrics'][independent]['MAE'][k]
            passed=baseline is not None and baseline>0 and candidate<=.9*baseline
            expected_comparisons.append(dict(baseline=primary,trait=trait,additive_MAE=candidate,baseline_MAE=baseline,
                absolute_MAE_difference=None if baseline is None else candidate-baseline,
                relative_gain=None if baseline is None or baseline==0 else (baseline-candidate)/baseline,passes_10_percent=passed))
            if not passed:reasons.append('MAE_hurdle_failed:'+primary+':'+trait)
        candidate=metrics['metrics']['additive']['misclassification'];baseline=metrics['metrics'][independent]['misclassification']
        if baseline is None or candidate>baseline:reasons.append('classification_error_worse_or_unavailable:'+primary)
    check.fields(score,{'MAE_comparisons':expected_comparisons,'failure_reasons':reasons},'predictive decisions')
    # Reconstruct all complete-square and full observed-component decisions.
    descriptions=evaluation['descriptions'];check.equal(len(descriptions),len(units),'description count')
    bygroup={};byseed={};strata={s:dict(selected_squares=0,**stats()) for s in ('isolated_observed_square','larger_observed_component')}
    geometry_cells=0;diagonals=0;scenarios=0
    for d,u,cov in zip(descriptions,units,coverage):
        check.fields(d,{'square_index':u['square_index'],'seed':u['seed'],'group':group[u['seed']],'isolated_observed_square':u['isolated_observed_square'],'eligible':cov['eligible']},'description')
        stratum='isolated_observed_square' if u['isolated_observed_square'] else 'larger_observed_component';strata[stratum]['selected_squares']+=1
        gs=bygroup.setdefault(group[u['seed']],stats());ss=byseed.setdefault(u['seed'],stats());tables=(gs,ss,strata[stratum])
        for table in tables:table['eligible_squares']+=int(cov['eligible'])
        if not cov['eligible']:
            check.fields(d,{'diagonals':[],'missing_columns':cov['missing_columns'],'unclassified_endpoint_IDs':u['diagonals']},'ineligible square');continue
        ids=u['coordinate_order_00_01_10_11'];component=cids[u['component_index']]
        local={j:point(rows[v]) for j,v in enumerate(ids)};full={v:point(rows[v]) for v in component}
        local_ranges={j:ranges(rows[v]) for j,v in enumerate(ids)};full_ranges={v:ranges(rows[v]) for v in component}
        cedges=[(a,b) for a,b in edges if a in full and b in full];delta=tolerances(full);idelta=range_tolerance(full_ranges)
        for k,trait in enumerate(('protein','growth')):
            inter=ref.interaction([local[j][k] for j in range(4)])
            for representation,e in (('raw',inter),('natural_log',inter['natural_log'])):
                a=d['additive_descriptions'][trait][representation]
                if e is None:check.equal(a,None,'excluded log')
                else:check.fields(a,{'contrast':e['contrast'],'projection':e['projection'],
                    'residual':[sign*e['contrast']/4 for sign in (1,-1,-1,1)]},'additive '+representation)
            check.equal(d['additive_descriptions'][trait]['log_exclusion_reason'],None if inter['natural_log'] is not None else 'nonpositive_value','log eligibility')
        check.equal(len(d['diagonals']),2,'diagonal count')
        for actual,ep in zip(d['diagonals'],ref.DIAGONALS):
            cep=tuple(ids[j] for j in ep);primary=classify(ref,local,ref.SQUARE_EDGES,ep)
            expected={'endpoint_IDs':list(cep),'primary':primary,'point_tolerance':delta,
                'local_alternative_tolerance':classify(ref,local,ref.SQUARE_EDGES,ep,delta),
                'observed_component':field_check(ref,full,cedges,cep),'observed_component_alternative_tolerance':field_check(ref,full,cedges,cep,delta),
                'local_range_certificate':bounds(ref,local_ranges,ref.SQUARE_EDGES,ep),'component_range_certificate':bounds(ref,full_ranges,cedges,cep),
                'interval_tolerance':idelta,'local_range_alternative_tolerance':bounds(ref,local_ranges,ref.SQUARE_EDGES,ep,idelta),
                'component_range_alternative_tolerance':bounds(ref,full_ranges,cedges,cep,idelta)}
            check.fields(actual,expected,'diagonal')
            geometry_cells+=check_geometry(check,ref,actual['threshold_geometry_tolerance_zero'],local,ep)
            tau=[min(local[j][k] for j in ep) for k in (0,1)]
            check.equal(actual['native_shortfalls'],{str(j):[max(0,tau[k]-local[j][k]) for k in (0,1)] for j in range(4) if j not in ep},'shortfalls')
            prep=[]
            for pcol in PROTEIN:
                for gcol in GROWTH:
                    lf={j:[rows[v][pcol]['value'],rows[v][gcol]['value']] for j,v in enumerate(ids)}
                    cf={v:[rows[v][pcol]['value'],rows[v][gcol]['value']] for v in component};dd=tolerances(cf)
                    prep.append({'protein_column':pcol,'growth_column':gcol,'local':classify(ref,lf,ref.SQUARE_EDGES,ep),
                        'component':field_check(ref,cf,cedges,cep),'tolerance':dd,'local_alternative_tolerance':classify(ref,lf,ref.SQUARE_EDGES,ep,dd),
                        'component_alternative_tolerance':field_check(ref,cf,cedges,cep,dd)})
            check.equal(actual['preparation_scenarios'],prep,'12 scenarios')
            consistent=all(p['local']['classification']=='joint_only_obstruction' for p in prep)
            check.equal(actual['preparation_consistent_joint_obstruction'],consistent,'scenario consistency')
            for table in tables:
                table['eligible_diagonals']+=1;table['joint_only_obstructions']+=int(primary['classification']=='joint_only_obstruction')
                table['preparation_consistent_diagonals']+=int(consistent)
                table['range_certified_diagonals']+=int(expected['local_range_certificate']['status']=='certified_joint_only_obstruction')
            diagonals+=1;scenarios+=12
    total={k:sum(g[k] for g in bygroup.values()) for k in stats()};total['eligible_related_seed_groups']=sum(g['eligible_squares']>0 for g in bygroup.values());fraction(total)
    for table in (bygroup,byseed,strata):
        for s in table.values():fraction(s)
    coverage_pass=total['eligible_squares']>=10 and total['eligible_related_seed_groups']>=8
    check.fields(evaluation,{'primary_summary':total,'summary_by_related_seed_group':bygroup,'summary_by_seed':byseed,
        'summary_by_observed_component_stratum':strata,'coverage_passed':coverage_pass,
        'arm_descriptive_milestone':coverage_pass and total['preparation_consistent_diagonals']>0,
        'failure_reasons':([] if coverage_pass else ['primary_coverage_below_10_squares_8_groups'])+([] if total['preparation_consistent_diagonals'] else ['no_preparation_consistent_joint_obstruction'])},'final summary')
    return {'status':'verified_complete_arm','primary_summary':total,'prediction_screen_passed':metrics['success'],
        'independent_prediction_metrics':metrics,'independent_summary_by_related_group':bygroup,
        'prediction_cases':metrics['n_cases'],'diagonals_verified':diagonals,'preparation_scenarios_verified':scenarios,
        'exact_critical_cells_verified':geometry_cells,'mean_checks':len(mean_checks)}

def run(authorization,root,out):
    out=Path(out)
    if out.exists():raise FileExistsError('Refuse existing output before authorization access')
    auth,paths,hashes,ref=authorize(Path(authorization),Path(root))
    out.mkdir(parents=True,exist_ok=False);check=Check()
    try:
        result=compare(paths,hashes,auth['arm'],ref,check)
        if any(sha(paths[r])!=h for r,h in hashes.items()):raise ValueError('Input changed during comparison')
        receipt=dict(result,arm=auth['arm'],data_scope=auth['data_scope'],input_hashes_before_after_equal=hashes,
            exact_leaf_comparisons=check.exact,numerical_identity_comparisons=check.numeric,
            numerical_identity_tolerance={'relative':TOLERANCE,'absolute':TOLERANCE},
            classifications_and_scientific_thresholds='exact; numerical comparison tolerance never changes scientific decisions',
            primary_implementation_imported=False,raw_archive_opened=False,
            limits=['Selected frozen exports only; does not authenticate raw processing, historical access order or biological uncertainty.',
                    'Supplemental amendment diagnostics require separate verification; this adapter covers the held core runner schema.'])
        write(out/'receipt.json',receipt);return receipt
    except Exception as exc:
        write(out/'failure.json',{'status':'mismatch_or_invalid_input','error':str(exc),'arm':auth['arm'],'data_scope':auth['data_scope'],
            'input_hashes':hashes,'exact_leaf_comparisons_before_stop':check.exact,'numerical_comparisons_before_stop':check.numeric})
        raise

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--authorization',required=True,type=Path)
    parser.add_argument('--root',required=True,type=Path);parser.add_argument('--out',required=True,type=Path)
    a=parser.parse_args();r=run(a.authorization,a.root,a.out);print(json.dumps({'status':r['status'],'arm':r['arm']}))

if __name__=='__main__':main()
