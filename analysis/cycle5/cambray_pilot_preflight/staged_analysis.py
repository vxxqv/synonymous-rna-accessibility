"""Portable, hash-gated Cambray pilot runner. Never run biological inputs unfrozen.

CLI: python staged_analysis.py predict|evaluate --authorization FILE --root DIR
     --arm development|reserved_validation --output NEW_DIRECTORY

Authorization schema Cambray_staged_v1 requires execution_authorized=true,
allowed_phases, and hash-pinned input roles protocol,split,mask,components,edges,
archive,graph_rules,masking,runner. Every input is {path,sha256}; paths resolve
against --root. The three code paths must equal the imported/local code files.
Use an absolute archive path to keep the 130 MB source outside the repository.
expected_records is 244000 for the authored source and may differ in toy tests.
Evaluation additionally pins prediction_seal. Reserved-arm execution additionally
pins development_decision and development_report. Decision JSON requires status
development_report_frozen, unchanged_validation_procedure=true, protocol_sha256,
report_sha256, and code_hashes keyed graph_rules,masking,runner.

Execution sequence: authorize development predict, freeze its seal, then authorize
development evaluate with that seal. Freeze the development report and an explicit
unchanged-procedure decision. Authorize reserved_validation predict with these
two inputs, freeze its seal, then authorize reserved_validation evaluate with that
seal and the same decision/report. Every invocation needs a fresh output directory.
The authorization is an external root decision, not a self-issued approval.

This module provides no network or Git access. Input metadata is not a security
sandbox: the separate root authorization and reviewed frozen inputs are required.
"""
if not __debug__:
    raise SystemExit('Disabled assertions are unsupported; refuse -O and -OO before input access')

import argparse
import csv
from fractions import Fraction
import hashlib
import io
import json
import math
from pathlib import Path
import re
import zipfile
import graph_rules as rules
import masking_v2 as masking

PROTEIN=['clean.lin.prot.repOne','clean.lin.prot.repTwo','clean.lin.prot.repThree','clean.lin.prot.repFour']
GROWTH=['growth.R60one','growth.R60two','growth.R60three']
REPLICATES=PROTEIN+GROWTH
AUTHOR_MEANS=['clean.lin.prot.mean','growth.mean.R60']
METHODS=['additive','adjacent_mean','three_corner_mean']
TRAITS=['protein','growth']
SQUARE=[[1,2],[0,3],[0,3],[1,2]]
BASE_ROLES=['protocol','split','mask','components','edges','archive','graph_rules','masking','runner']


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def write_json(path,value):
    def encode(x):
        if isinstance(x,Fraction):return {'numerator':x.numerator,'denominator':x.denominator}
        raise TypeError(type(x).__name__)
    with Path(path).open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,sort_keys=True,allow_nan=False,default=encode);stream.write('\n')


def authorize(path,root,phase,arm):
    auth=json.loads(Path(path).read_text(encoding='utf-8'))
    if auth.get('schema')!='Cambray_staged_v1' or auth.get('execution_authorized') is not True or phase not in auth.get('allowed_phases',[]):
        raise ValueError('Execution authorization absent or phase refused')
    if arm not in auth.get('allowed_arms',[]):raise ValueError('Arm not authorized')
    roles=BASE_ROLES+(['prediction_seal'] if phase=='evaluate' else [])
    if arm=='reserved_validation':roles+=['development_decision','development_report']
    paths={};hashes={}
    for role in roles:
        entry=auth['inputs'][role];p=Path(entry['path'])
        if not p.is_absolute():p=Path(root)/p
        p=p.resolve()
        actual=sha(p)
        if actual!=entry['sha256']:raise ValueError('Input hash mismatch: '+role)
        paths[role]=p;hashes[role]=actual
    for role,module in [('runner',Path(__file__)),('graph_rules',Path(rules.__file__)),('masking',Path(masking.__file__))]:
        if paths[role]!=module.resolve():raise ValueError('Pinned code is not executing code: '+role)
    if arm=='reserved_validation':
        decision=json.loads(paths['development_decision'].read_text())
        if decision.get('status')!='development_report_frozen' or decision.get('unchanged_validation_procedure') is not True:
            raise ValueError('Development decision does not authorize unchanged validation procedure')
        if decision.get('report_sha256')!=hashes['development_report'] or decision.get('protocol_sha256')!=hashes['protocol']:
            raise ValueError('Development report/protocol hash mismatch')
        if decision.get('code_hashes')!={r:hashes[r] for r in ('graph_rules','masking','runner')}:
            raise ValueError('Code changed since development decision')
    return auth,paths,hashes


def row_identifier(line):
    token=line.partition(b',')[0]
    if re.fullmatch(rb'"[A-Za-z0-9_.-]+"',token) is None:raise ValueError('Unexpected row ID syntax')
    return token[1:-1].decode('ascii')


def selected_fields(line,positions):
    """Tokenize bytes only as far as requested fields; decode no other values.

    The authenticated author format has one physical line per row. Quoted commas
    and doubled quotes are handled. Selected multiline fields are unsupported.
    """
    wanted=set(positions);last=max(wanted);column=0;start=0;quoted=False;i=0;found={}
    while i<=len(line):
        current=line[i] if i<len(line) else None
        if current==34:
            if quoted and i+1<len(line) and line[i+1]==34:i+=2;continue
            quoted=not quoted
        if (current==44 and not quoted) or current is None or (current in (10,13) and not quoted):
            if column in wanted:found[column]=line[start:i]
            if column==last:return found
            column+=1;start=i+1
            if current is None or current in (10,13):break
        i+=1
    raise ValueError('Selected field missing or malformed CSV boundary')


def parse_numeric(token):
    if token.startswith(b'"') and token.endswith(b'"'):token=token[1:-1].replace(b'""',b'"')
    text=token.decode('ascii').strip()
    if text in ('','NA','NaN','nan','NULL'):return {'value':None,'missing_reason':text or 'empty','raw':text}
    try:value=float(text)
    except ValueError as exc:raise ValueError('Invalid requested numeric field') from exc
    if not math.isfinite(value):return {'value':None,'missing_reason':'nonfinite','raw':text}
    return {'value':value,'missing_reason':None,'raw':text}


def selective_reader(archive,allowed_ids,fields,expected_records=None):
    """Identity first. Hidden rows' phenotype bytes are never tokenized/decoded."""
    allowed=set(allowed_ids);rows={};seen=set();total=0
    if not set(fields)<=set(REPLICATES+AUTHOR_MEANS):raise ValueError('Undeclared phenotype field requested')
    with zipfile.ZipFile(archive) as zf:
        candidates=[n for n in zf.namelist() if n.endswith('.csv') and not n.startswith('__MACOSX/')]
        if len(candidates)!=1:raise ValueError('Ambiguous CSV member')
        with zf.open(candidates[0]) as stream:
            header=next(csv.reader(io.StringIO(next(stream).decode('utf-8'))))
            if len(set(header))!=len(header):raise ValueError('Duplicate column names')
            positions={name:header.index(name) for name in fields}
            for line in stream:
                identifier=row_identifier(line);total+=1
                if identifier in seen:raise ValueError('Duplicate source ID')
                seen.add(identifier)
                if identifier not in allowed:continue
                selected=selected_fields(line,positions.values())
                rows[identifier]={name:parse_numeric(selected[pos]) for name,pos in positions.items()}
    if expected_records is not None and total!=expected_records:raise ValueError('Source row count mismatch')
    if set(rows)!=allowed:raise ValueError('Selected identity absent from source')
    return rows,{'source_records':total,'selected_records':len(rows),'decoded_fields':list(fields),'other_row_phenotype_values_decoded':0}


def load_metadata(paths,arm):
    split=json.loads(paths['split'].read_text());mask=json.loads(paths['mask'].read_text())
    units=[u for u in split['units'] if u['arm']==arm]
    fixed={u['withheld_ID'] for u in units}
    if fixed!=set(mask['arms'][arm]['masked_IDs']):raise ValueError('Union mask differs from every original hashed target')
    mask_units={u['square_index']:u for u in mask['units'] if u['arm']==arm}
    if set(mask_units)!={u['square_index'] for u in units}:raise ValueError('Mask/unit discrepancy')
    for unit in units:
        m=mask_units[unit['square_index']]
        eligible=len(set(unit['coordinate_order_00_01_10_11'])&fixed)==1
        if m['prediction_eligible_from_metadata']!=eligible or m['coordinate_order_00_01_10_11']!=unit['coordinate_order_00_01_10_11']:
            raise ValueError('Mask eligibility changed')
    groups={seed:','.join(group) for group in split['seed_groups'] for seed in group}
    components=[json.loads(line) for line in paths['components'].read_text().splitlines()]
    return split,units,mask_units,fixed,groups,components


def complete_mean(row,columns):
    values=[row[c]['value'] for c in columns]
    if any(v is None for v in values):return None
    result=math.fsum(values)/len(values)
    if not math.isfinite(result):raise ValueError('Nonfinite derived mean')
    return result


def point(row):return [complete_mean(row,columns) for columns in (PROTEIN,GROWTH)]


def ranges(row):
    result=[]
    for columns in (PROTEIN,GROWTH):
        values=[row[c]['value'] for c in columns]
        result.append(None if any(v is None for v in values) else [min(values),max(values)])
    return result


def missing_columns(rows,ids):
    return [{'ID':i,'column':c,'reason':rows[i][c]['missing_reason']} for i in ids for c in REPLICATES if rows[i][c]['value'] is None]


def predict_stage(paths,hashes,auth,arm,output):
    split,units,mask_units,mask,groups,components=load_metadata(paths,arm)
    ids={i for u in units if mask_units[u['square_index']]['prediction_eligible_from_metadata'] for i in u['coordinate_order_00_01_10_11'] if i not in mask}
    if ids&mask:raise AssertionError('Masked ID entered visible extraction')
    rows,reader_receipt=selective_reader(paths['archive'],ids,REPLICATES,auth.get('expected_records'))
    predictions=[]
    for unit in units:
        result={'square_index':unit['square_index'],'seed':unit['seed'],'group':groups[unit['seed']],'withheld_ID':unit['withheld_ID']}
        if not mask_units[unit['square_index']]['prediction_eligible_from_metadata']:
            result.update(status='excluded_shared_corner_mask',reason='multiple_globally_masked_corners');predictions.append(result);continue
        coordinate=unit['coordinate_order_00_01_10_11'];visible=[i for i in coordinate if i not in mask]
        missing=missing_columns(rows,visible)
        if missing:
            result.update(status='missing_known_replicates',missing_columns=missing);predictions.append(result);continue
        known={i:point(rows[i]) for i in visible}
        per_trait=[masking.predict_from_visible_ids(coordinate,mask,{i:known[i][k] for i in visible}) for k in range(2)]
        predictions_by_method={name:[per_trait[k]['predictions'][name] for k in range(2)] for name in METHODS}
        missing_index=coordinate.index(unit['withheld_ID']);endpoints=(missing_index^1,missing_index^2)
        classifications={}
        for name in METHODS:
            field=[predictions_by_method[name] if j==missing_index else known[i] for j,i in enumerate(coordinate)]
            label=rules.classify(SQUARE,field,endpoints)['classification']
            classifications[name]={'classification':label,'joint_only_obstruction':label=='joint_only_obstruction'}
        result.update(status='predicted',visible_IDs=visible,known_corner_means=known,predictions=predictions_by_method,
                      classifications=classifications,classification_endpoint_IDs=[coordinate[j] for j in endpoints])
        predictions.append(result)
    write_json(output/'known_inputs.json',{'arm':arm,'rows':rows,'reader':reader_receipt})
    write_json(output/'predictions.json',{'arm':arm,'records':predictions,'masked_IDs':sorted(mask),'target_finiteness_examined':False,
              'all_source_and_metadata_hashes':{r:hashes[r] for r in BASE_ROLES}})
    seal={'schema':'Cambray_prediction_seal_v1','arm':arm,'status':'predictions_sealed_before_target_access',
          'input_hashes':{r:hashes[r] for r in BASE_ROLES},'files':{name:sha(output/name) for name in ('known_inputs.json','predictions.json')},
          'target_finiteness_examined':False,'author_means_examined':False,'component_outcomes_examined':False}
    write_json(output/'seal.json',seal)
    return {'status':seal['status'],'arm':arm,'prediction_records':len(predictions),'seal_sha256':sha(output/'seal.json')}


def verify_prediction_seal(paths,hashes,arm):
    path=paths['prediction_seal'];seal=json.loads(path.read_text())
    if seal.get('schema')!='Cambray_prediction_seal_v1' or seal.get('arm')!=arm or seal.get('status')!='predictions_sealed_before_target_access':
        raise ValueError('Invalid prediction seal')
    if seal.get('input_hashes')!={r:hashes[r] for r in BASE_ROLES}:raise ValueError('Inputs/source changed after predictions')
    if any(seal.get(k) is not False for k in ('target_finiteness_examined','author_means_examined','component_outcomes_examined')):
        raise ValueError('Prediction stage reports premature exposure')
    if set(seal['files'])!={'known_inputs.json','predictions.json'}:raise ValueError('Unexpected sealed files')
    for name,digest in seal['files'].items():
        if sha(path.parent/name)!=digest:raise ValueError('Sealed numerical payload changed')
    return json.loads((path.parent/'predictions.json').read_text()),json.loads((path.parent/'known_inputs.json').read_text())


def graph_for(ids,adjacency):
    index={identifier:i for i,identifier in enumerate(ids)}
    return [sorted(index[v] for v in adjacency[identifier] if v in index) for identifier in ids]


def field_check(graph,field,endpoints,tolerance):
    return {'complete':all(v is not None for row in field for v in row),
            'point':rules.classify(graph,field,endpoints,tolerance=tolerance) if all(v is not None for row in field for v in row) else None,
            'bounds':rules.missing_vertex_bounds(graph,field,endpoints,tolerance=tolerance)}


def interval_tolerance(intervals):
    values=[]
    for row in intervals:
        values.append([None if interval is None else max(abs(interval[0]),abs(interval[1])) for interval in row])
    return rules.component_tolerances(values)


def assess_diagonal(unit,rows,component_ids,component_graph,endpoints):
    ids=unit['coordinate_order_00_01_10_11'];c_end=tuple(component_ids.index(ids[i]) for i in endpoints)
    local=[point(rows[i]) for i in ids];full=[point(rows[i]) for i in component_ids]
    tol=rules.component_tolerances(full)
    local_intervals=[ranges(rows[i]) for i in ids];full_intervals=[ranges(rows[i]) for i in component_ids]
    itol=interval_tolerance(full_intervals)
    result={'endpoint_IDs':[ids[i] for i in endpoints],'primary':rules.classify(SQUARE,local,endpoints),
            'point_tolerance':tol,'local_alternative_tolerance':rules.classify(SQUARE,local,endpoints,tolerance=tol),
            'observed_component':field_check(component_graph,full,c_end,0),
            'observed_component_alternative_tolerance':field_check(component_graph,full,c_end,tol),
            'threshold_geometry_tolerance_zero':rules.threshold_geometry(SQUARE,local,endpoints),
            'native_shortfalls':rules.native_shortfalls(local,endpoints),
            'local_range_certificate':rules.interval_certificate(SQUARE,local_intervals,endpoints),
            'component_range_certificate':rules.interval_certificate(component_graph,full_intervals,c_end),
            'interval_tolerance':itol,
            'local_range_alternative_tolerance':rules.interval_certificate(SQUARE,local_intervals,endpoints,tolerance=itol),
            'component_range_alternative_tolerance':rules.interval_certificate(component_graph,full_intervals,c_end,tolerance=itol),
            'preparation_scenarios':[]}
    for protein in PROTEIN:
        for growth in GROWTH:
            field=lambda identifiers:[[rows[i][protein]['value'],rows[i][growth]['value']] for i in identifiers]
            lfield=field(ids);cfield=field(component_ids);delta=rules.component_tolerances(cfield)
            result['preparation_scenarios'].append({'protein_column':protein,'growth_column':growth,
                'local':rules.classify(SQUARE,lfield,endpoints),'component':field_check(component_graph,cfield,c_end,0),
                'tolerance':delta,'local_alternative_tolerance':rules.classify(SQUARE,lfield,endpoints,tolerance=delta),
                'component_alternative_tolerance':field_check(component_graph,cfield,c_end,delta)})
    result['preparation_consistent_joint_obstruction']=all(s['local']['classification']=='joint_only_obstruction' for s in result['preparation_scenarios'])
    return result


def grouped_mean(records,value):
    grouped={}
    for r in records:grouped.setdefault(r['group'],[]).append(value(r))
    per_group={g:math.fsum(v)/len(v) for g,v in grouped.items()}
    return {'value':math.fsum(per_group.values())/len(per_group) if per_group else None,'groups':per_group}


def score_predictions(predictions,unit_map,rows,mask_units):
    records=[]
    for prediction in predictions['records']:
        unit=unit_map[prediction['square_index']];ids=unit['coordinate_order_00_01_10_11']
        record={'square_index':unit['square_index'],'seed':unit['seed'],'group':prediction['group'],'status':prediction['status']}
        if prediction['status']!='predicted':
            record['reason']=prediction.get('reason',prediction.get('missing_columns'));records.append(record);continue
        missing=missing_columns(rows,ids)
        if missing:
            record.update(status='unscorable_missing_target_or_square',missing_columns=missing);records.append(record);continue
        observed=point(rows[prediction['withheld_ID']]);position=ids.index(prediction['withheld_ID']);endpoints=(position^1,position^2)
        true_class=rules.classify(SQUARE,[point(rows[i]) for i in ids],endpoints)['classification']=='joint_only_obstruction'
        record.update(status='scored',observed_target=observed,true_joint_only_obstruction=true_class,methods={})
        for method in METHODS:
            estimate=prediction['predictions'][method];pred_class=prediction['classifications'][method]['joint_only_obstruction']
            error=[p-y for p,y in zip(estimate,observed)]
            record['methods'][method]={'predicted_target':estimate,'signed_error':error,'absolute_error':list(map(abs,error)),
                'squared_error':[v*v for v in error],'predicted_joint_only_obstruction':pred_class,'classification_error':int(pred_class!=true_class)}
        records.append(record)
    scored=[r for r in records if r['status']=='scored'];summary={}
    for method in METHODS:
        continuous={}
        for k,trait in enumerate(TRAITS):
            mae=grouped_mean(scored,lambda r:r['methods'][method]['absolute_error'][k])
            mse=grouped_mean(scored,lambda r:r['methods'][method]['squared_error'][k])
            continuous[trait]={'MAE':mae['value'],'RMSE':math.sqrt(mse['value']) if mse['value'] is not None else None,'group_MAE':mae['groups'],'group_MSE':mse['groups']}
        confusion={label:0 for label in ('TP','FN','TN','FP')}
        for r in scored:
            true=r['true_joint_only_obstruction'];pred=r['methods'][method]['predicted_joint_only_obstruction']
            confusion[('TP' if pred else 'FN') if true else ('FP' if pred else 'TN')]+=1
        positives=confusion['TP']+confusion['FN'];negatives=confusion['TN']+confusion['FP']
        balanced=(confusion['TP']/positives+confusion['TN']/negatives)/2 if positives and negatives else None
        summary[method]={'continuous':continuous,'pooled_confusion':confusion,'pooled_balanced_accuracy':balanced,
            'group_balanced_classification_error':grouped_mean(scored,lambda r:r['methods'][method]['classification_error'])}
    groups={r['group'] for r in scored};both_classes=len({r['true_joint_only_obstruction'] for r in scored})==2
    comparisons=[];reasons=[]
    if len(scored)<10 or len(groups)<8:reasons.append('prediction_coverage_below_10_squares_8_groups')
    if not both_classes:reasons.append('both_true_classes_not_observed')
    for baseline in METHODS[1:]:
        for trait in TRAITS:
            candidate=summary['additive']['continuous'][trait]['MAE'];base=summary[baseline]['continuous'][trait]['MAE']
            passed=base is not None and base>0 and candidate<=.9*base
            comparisons.append({'baseline':baseline,'trait':trait,'additive_MAE':candidate,'baseline_MAE':base,
                'absolute_MAE_difference':None if base is None else candidate-base,
                'relative_gain':None if base is None or base==0 else (base-candidate)/base,'passes_10_percent':passed})
            if not passed:reasons.append('MAE_hurdle_failed:'+baseline+':'+trait)
        candidate=summary['additive']['group_balanced_classification_error']['value'];base=summary[baseline]['group_balanced_classification_error']['value']
        if base is None or candidate>base:reasons.append('classification_error_worse_or_unavailable:'+baseline)
    return {'records':records,'scored_squares':len(scored),'related_seed_groups':len(groups),'both_true_classes':both_classes,
        'summaries':summary,'MAE_comparisons':comparisons,'discrimination_screen_passed':not reasons,'failure_reasons':reasons,
        'scope':'Conditional missing-corner prediction. Shared targets and preparation dependence retained; no inference intervals.'}


def evaluate_stage(paths,hashes,auth,arm,output):
    #All seal checks precede calling the phenotype reader or target-finiteness test.
    predictions,known=verify_prediction_seal(paths,hashes,arm)
    split,units,mask_units,mask,groups,components=load_metadata(paths,arm)
    component_ids={u['component_index']:sorted(components[u['component_index']]['IDs']) for u in units}
    all_ids={v for ids in component_ids.values() for v in ids}
    rows,reader=selective_reader(paths['archive'],all_ids,REPLICATES+AUTHOR_MEANS,auth.get('expected_records'))
    for identifier,row in known['rows'].items():
        if {c:rows[identifier][c] for c in REPLICATES}!=row:raise ValueError('Previously visible input changed')
    mean_checks=[];holds=[]
    for identifier in sorted(rows):
        for k,columns in enumerate((PROTEIN,GROWTH)):
            value=complete_mean(rows[identifier],columns);archived=rows[identifier][AUTHOR_MEANS[k]]['value']
            if value is None:status='not_comparable_missing_replicates'
            elif archived is None:status='missing_author_mean';holds.append((identifier,TRAITS[k],status))
            elif math.isclose(value,archived,rel_tol=1e-8,abs_tol=1e-10):status='compatible'
            else:status='mean_discrepancy';holds.append((identifier,TRAITS[k],status))
            mean_checks.append({'ID':identifier,'trait':TRAITS[k],'computed':value,'author':archived,'difference':None if value is None or archived is None else value-archived,
                'tolerance':None if value is None or archived is None else max(1e-10,1e-8*max(abs(value),abs(archived))),'status':status})
    coverage=[]
    for unit in units:
        missing=missing_columns(rows,unit['coordinate_order_00_01_10_11'])
        coverage.append({'square_index':unit['square_index'],'seed':unit['seed'],'group':groups[unit['seed']],
            'eligible':not missing,'missing_columns':missing,'all_diagonals_retained':unit['diagonals'],
            'prediction_metadata_eligible':mask_units[unit['square_index']]['prediction_eligible_from_metadata']})
    write_json(output/'unsealed_inputs.json',{'arm':arm,'rows':rows,'reader':reader})
    base={'arm':arm,'mean_diagnostics':mean_checks,'coverage':coverage,'source_hashes':hashes,'prediction_seal_sha256':hashes['prediction_seal']}
    if holds:
        base.update(status='held_for_author_mean_provenance',holds=holds,interpretation_permitted=False)
        write_json(output/'evaluation.json',base);return {'status':base['status'],'arm':arm,'holds':len(holds)}
    adjacency={v:set() for v in all_ids}
    with paths['edges'].open(encoding='utf-8') as stream:
        for edge in csv.DictReader(stream,delimiter='\t'):
            a,b=edge['left_ID'],edge['right_ID']
            if a in all_ids and b in all_ids:adjacency[a].add(b);adjacency[b].add(a)
    descriptions=[]
    for unit,cov in zip(units,coverage):
        result={'square_index':unit['square_index'],'seed':unit['seed'],'group':groups[unit['seed']],
                'isolated_observed_square':unit['isolated_observed_square'],'eligible':cov['eligible'],'diagonals':[]}
        if not cov['eligible']:
            result['missing_columns']=cov['missing_columns'];result['unclassified_endpoint_IDs']=unit['diagonals'];descriptions.append(result);continue
        cid=unit['component_index'];ids=component_ids[cid];graph=graph_for(ids,adjacency)
        field=[point(rows[v]) for v in unit['coordinate_order_00_01_10_11']]
        result['additive_descriptions']={}
        for k,trait in enumerate(TRAITS):
            values=[row[k] for row in field]
            result['additive_descriptions'][trait]={'raw':rules.additive_projection(values),
                'natural_log':rules.additive_projection([math.log(v) for v in values]) if all(v>0 for v in values) else None,
                'log_exclusion_reason':None if all(v>0 for v in values) else 'nonpositive_value'}
        for endpoints in ((0,3),(1,2)):result['diagonals'].append(assess_diagonal(unit,rows,ids,graph,endpoints))
        descriptions.append(result)
    by_group={};by_seed={};by_stratum={name:{'selected_squares':0,'eligible_squares':0,'eligible_diagonals':0,
        'joint_only_obstructions':0,'preparation_consistent_diagonals':0,'range_certified_diagonals':0}
        for name in ('isolated_observed_square','larger_observed_component')}
    for result in descriptions:
        stratum='isolated_observed_square' if result['isolated_observed_square'] else 'larger_observed_component'
        by_stratum[stratum]['selected_squares']+=1
        for table,key in ((by_group,result['group']),(by_seed,result['seed']),(by_stratum,stratum)):
            stats=table.setdefault(key,{'eligible_squares':0,'eligible_diagonals':0,'joint_only_obstructions':0,'preparation_consistent_diagonals':0,'range_certified_diagonals':0})
            stats['eligible_squares']+=int(result['eligible'])
            for diagonal in result['diagonals']:
                stats['eligible_diagonals']+=1
                stats['joint_only_obstructions']+=diagonal['primary']['classification']=='joint_only_obstruction'
                stats['preparation_consistent_diagonals']+=diagonal['preparation_consistent_joint_obstruction']
                stats['range_certified_diagonals']+=diagonal['local_range_certificate']['status']=='certified_joint_only_obstruction'
    totals={key:sum(g[key] for g in by_group.values()) for key in ('eligible_squares','eligible_diagonals','joint_only_obstructions','preparation_consistent_diagonals','range_certified_diagonals')}
    totals['eligible_related_seed_groups']=sum(g['eligible_squares']>0 for g in by_group.values())
    totals['joint_only_fraction']=totals['joint_only_obstructions']/totals['eligible_diagonals'] if totals['eligible_diagonals'] else None
    for table in (by_group,by_seed,by_stratum):
        for stats in table.values():stats['joint_only_fraction']=stats['joint_only_obstructions']/stats['eligible_diagonals'] if stats['eligible_diagonals'] else None
    coverage_pass=totals['eligible_squares']>=10 and totals['eligible_related_seed_groups']>=8
    scores=score_predictions(predictions,{u['square_index']:u for u in units},rows,mask_units)
    base.update(status='arm_evaluation_complete',interpretation_permitted=True,descriptions=descriptions,primary_summary=totals,
        summary_by_related_seed_group=by_group)
    #The two-stage parser permits no validation adaptation; interpretation is
    #limited to this arm until the root combines the separately sealed decisions.
    base['summary_by_seed']=by_seed;base['coverage_passed']=coverage_pass
    base['summary_by_observed_component_stratum']=by_stratum
    base['arm_descriptive_milestone']=coverage_pass and totals['preparation_consistent_diagonals']>0
    base['prediction_evaluation']=scores
    base['failure_reasons']=([] if coverage_pass else ['primary_coverage_below_10_squares_8_groups'])+([] if totals['preparation_consistent_diagonals'] else ['no_preparation_consistent_joint_obstruction'])
    write_json(output/'evaluation.json',base)
    write_json(output/'evaluation_receipt.json',{'status':base['status'],'arm':arm,'files':{n:sha(output/n) for n in ('unsealed_inputs.json','evaluation.json')},'input_hashes':hashes})
    return {'status':base['status'],'arm':arm,'primary_summary':totals,'prediction_screen_passed':scores['discrimination_screen_passed']}


def run(phase,authorization,root,arm,output):
    if phase not in ('predict','evaluate') or arm not in ('development','reserved_validation'):
        raise ValueError('Unsupported stage or arm')
    output=Path(output)
    if output.exists():raise FileExistsError('Refuse existing output directory')
    auth,paths,hashes=authorize(authorization,root,phase,arm)
    output.mkdir(parents=True,exist_ok=False)
    return predict_stage(paths,hashes,auth,arm,output) if phase=='predict' else evaluate_stage(paths,hashes,auth,arm,output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase',choices=['predict','evaluate']);parser.add_argument('--authorization',required=True)
    parser.add_argument('--root',required=True);parser.add_argument('--arm',required=True,choices=['development','reserved_validation'])
    parser.add_argument('--output',required=True);args=parser.parse_args()
    print(json.dumps(run(args.phase,args.authorization,args.root,args.arm,args.output),indent=2,allow_nan=False))

if __name__=='__main__':main()
