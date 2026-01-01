"""Independent DEVELOPMENT-only supplement arithmetic, without author imports.

Consumes only the pinned evaluation, diagnostics and sequence-only split.
Stored predictions/targets define exact binary rational inputs to the loss oracle.
Underlying graph classifications are audited separately by another review lane.
"""
if not __debug__:
    raise RuntimeError('Optimization disabled checks are refused')

from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import platform

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
REPO = WORKSPACE / 'outputs/synonymous-rna-accessibility'
PINS = {
 'evaluation': ('results/cycle5/development/evaluation/evaluation.json', '756ff513c7c4f628110b07de19c0608ef5b1854da548b0227bca44db75e8f8fa'),
 'diagnostics': ('results/cycle5/development/diagnostics/diagnostics.json', 'c5a311f6470481ca6c0de474d7e80fdd348e100334822def01e2a342655f66db'),
 'split': ('analysis/cycle5/cambray_pilot_preflight/split_manifest.json', '394a0c1aa23785c290b31c147d199f406cf46b269daa9159a98c419398b05a4e'),
 'protocol': ('analysis/cycle5/cambray_pilot/protocol_v1.md', '0f51ccc808559ea69df3bdf294520d81babb67bf14ca92ce2abb6fcc665760ac'),
}
SCRIPT_SHA = '5e8c1a5638d0d2746f89f12cf26d209f11930f490df99c16a79eda3cf2628850'
METHODS = ('additive', 'adjacent_mean', 'three_corner_mean')
TRAITS = ('protein', 'growth')
METRICS = ('joint_only', 'preparation_consistent', 'range_certificate')
CLASS = ('jointly_connected', 'joint_only_obstruction', 'single_trait_disconnection', 'ineligible_endpoints', 'point_unavailable')
CERT = ('certified_jointly_connected', 'certified_joint_only_obstruction', 'certified_single_trait_disconnection', 'inconclusive', 'ineligible_definite_endpoints', 'certificate_unavailable')
BYPASS = ('certified_observed_component_bypass', 'observed_component_optimistic_disconnection', 'observed_component_bypass_unresolved', 'not_a_local_joint_only_obstruction', 'local_classification_unavailable')
CHECKS = Counter()
MAX_ABS = 0.0


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, label):
    CHECKS['conditions'] += 1
    if not condition:
        raise RuntimeError(label)


def equal(got, want, label):
    global MAX_ABS
    if isinstance(want, dict):
        require(isinstance(got, dict) and set(got) == set(want), label + ': keys')
        for key in want:
            equal(got[key], want[key], label + '/' + str(key))
    elif isinstance(want, list):
        require(isinstance(got, list) and len(got) == len(want), label + ': list size')
        for i, value in enumerate(want):
            equal(got[i], value, label + '/' + str(i))
    elif isinstance(want, bool) or want is None or isinstance(want, str):
        CHECKS['categorical_leaves'] += 1
        require(type(got) is type(want) and got == want, label + ': categorical')
    elif isinstance(want, int):
        CHECKS['integer_leaves'] += 1
        require(type(got) is int and got == want, label + ': integer')
    else:
        CHECKS['numeric_leaves'] += 1
        wanted = float(want)
        require(isinstance(got, (int, float)) and not isinstance(got, bool), label + ': numeric type')
        require(math.isfinite(got) and math.isfinite(wanted), label + ': nonfinite')
        MAX_ABS = max(MAX_ABS, abs(got-wanted))
        require(math.isclose(got, wanted, rel_tol=1e-12, abs_tol=1e-12), label + ': numerical mismatch')


def meanq(values):
    return sum(values, Q(0))/len(values) if values else None


def descriptive(rows, groups):
    group_counts = {}
    for group in groups:
        squares = [r for r in rows if r['group'] == group]
        diagonals = [d for r in squares for d in r['diagonals']]
        n = len(diagonals)
        counts = {
            'eligible_squares': sum(int(r['eligible']) for r in squares),
            'eligible_diagonals': n,
            'joint_only': sum(int(d['primary']['classification'] == 'joint_only_obstruction') for d in diagonals),
            'preparation_consistent': sum(int(all(s['local']['classification'] == 'joint_only_obstruction' for s in d['preparation_scenarios'])) for d in diagonals),
            'range_certificate': sum(int(d['local_range_certificate']['status'] == 'certified_joint_only_obstruction') for d in diagonals),
        }
        counts['fractions'] = {m: Q(counts[m], n) if n else None for m in METRICS}
        group_counts[group] = counts
    keys = ('eligible_squares', 'eligible_diagonals') + METRICS
    total = {key: sum(v[key] for v in group_counts.values()) for key in keys}
    total['fractions'] = {m: Q(total[m], total['eligible_diagonals']) if total['eligible_diagonals'] else None for m in METRICS}
    active = [g for g in groups if group_counts[g]['eligible_diagonals']]
    coverage = total['eligible_squares'] >= 10 and len(active) >= 8
    return {'raw_local': total, 'group_denominators': group_counts,
        'selected_groups': len(groups), 'contributing_groups': len(active),
        'zero_eligible_groups': [g for g in groups if g not in active],
        'equal_contributing_group_mean_fractions': {m: meanq([group_counts[g]['fractions'][m] for g in active]) for m in METRICS},
        'coverage_hurdle_passed': coverage,
        'descriptive_existence_milestone': coverage and total['preparation_consistent'] > 0}


def prediction(rows, groups, targets):
    selected = [r for r in rows if r['group'] in groups]
    active = sorted({r['group'] for r in selected})
    # Compute explicit rational case weights. No saved aggregate or saved loss is used.
    def loss(method, k, kind, group=None):
        use = [r for r in selected if group is None or r['group'] == group]
        if not use:
            return None
        result = Q(0)
        for row in use:
            same = [r for r in use if r['group'] == row['group']]
            if targets:
                count = len({r['withheld_ID'] for r in same}) * sum(r['withheld_ID'] == row['withheld_ID'] for r in same)
            else:
                count = len(same)
            if group is None:
                count *= len(active)
            if kind == 'classification':
                error = Q(int(row['methods'][method]['predicted_joint_only_obstruction'] != row['true_joint_only_obstruction']))
            else:
                error = Q(row['methods'][method]['predicted_target'][k])-Q(row['observed_target'][k])
                error = abs(error) if kind == 'absolute' else error*error
            result += error/count
        return result
    summaries = {}
    for method in METHODS:
        confusion = dict.fromkeys(('TP', 'TN', 'FP', 'FN'), 0)
        for r in selected:
            truth, pred = r['true_joint_only_obstruction'], r['methods'][method]['predicted_joint_only_obstruction']
            confusion[('TP' if truth else 'FP') if pred else ('FN' if truth else 'TN')] += 1
        positive, negative = confusion['TP']+confusion['FN'], confusion['TN']+confusion['FP']
        continuous = {}
        for k, trait in enumerate(TRAITS):
            mse = loss(method, k, 'square')
            continuous[trait] = {'MAE': loss(method, k, 'absolute'),
                'RMSE': math.sqrt(float(mse)) if mse is not None else None,
                'group_MAE': {g: loss(method, k, 'absolute', g) for g in groups},
                'group_MSE': {g: loss(method, k, 'square', g) for g in groups}}
        summaries[method] = {'continuous': continuous,
            'group_balanced_classification_error': loss(method, 0, 'classification'),
            'group_classification_error': {g: loss(method, 0, 'classification', g) for g in groups},
            'pooled_case_confusion': confusion,
            'pooled_case_balanced_accuracy': (Q(confusion['TP'],positive)+Q(confusion['TN'],negative))/2 if positive and negative else None}
    both = len({r['true_joint_only_obstruction'] for r in selected}) == 2
    reasons = []
    if len(selected) < 10 or len(active) < 8:
        reasons.append('prediction_coverage_below_10_squares_8_groups')
    if not both:
        reasons.append('both_true_classes_not_observed')
    comparisons = []
    for baseline in METHODS[1:]:
        for trait in TRAITS:
            a = summaries['additive']['continuous'][trait]['MAE']
            b = summaries[baseline]['continuous'][trait]['MAE']
            passed = b is not None and b > 0 and a <= Q(9,10)*b
            comparisons.append({'baseline': baseline, 'trait': trait, 'additive_MAE': a, 'baseline_MAE': b,
                'absolute_MAE_difference': a-b if b is not None else None,
                'relative_gain': (b-a)/b if b is not None and b > 0 else None,
                'passes_10_percent': passed})
            if not passed:
                reasons.append('MAE_hurdle_failed:'+baseline+':'+trait)
        a = summaries['additive']['group_balanced_classification_error']
        b = summaries[baseline]['group_balanced_classification_error']
        if b is None or a > b:
            reasons.append('classification_error_worse_or_unavailable:'+baseline)
    return {'weighting': 'group_target_case_loss' if targets else 'group_case_loss',
        'scored_cases': len(selected), 'selected_groups': len(groups), 'contributing_groups': len(active),
        'distinct_targets': len({r['withheld_ID'] for r in selected}), 'both_true_classes': both,
        'groups': {g: {'cases':sum(r['group']==g for r in selected), 'distinct_targets':len({r['withheld_ID'] for r in selected if r['group']==g})} for g in groups},
        'summaries': summaries, 'MAE_comparisons': comparisons,
        'discrimination_screen_passed': not reasons, 'failure_reasons': reasons}


def flattened(value):
    omitted = {'groups','group_denominators','group_MAE','group_MSE','group_classification_error','omitted_group'}
    out = {}; stack = [('', value)]
    while stack:
        name, item = stack.pop()
        if isinstance(item, dict):
            stack.extend(((name+'/'+k).lstrip('/'),v) for k,v in item.items() if k not in omitted)
        elif isinstance(item, list):
            stack.extend(((name+'/'+str(i)).lstrip('/'),v) for i,v in enumerate(item))
        elif item is None or isinstance(item,(int,float,bool,Q)):
            out[name] = item
    return out


def classification(field):
    return field['classification'] if field is not None else 'point_unavailable'


def certificate(field):
    return field['status'] if field is not None else 'certificate_unavailable'


def bypass(local, comp):
    if local is None:
        return BYPASS[4]
    if local['classification'] != 'joint_only_obstruction':
        return BYPASS[3]
    b = comp.get('bounds') if comp else None
    if b and b.get('joint_pessimistic') is True:
        return BYPASS[0]
    if b and b.get('joint_optimistic') is False:
        return BYPASS[1]
    return BYPASS[2]


def counts(values, categories):
    c = Counter(values)
    require(set(c) <= set(categories), 'Unexpected component category')
    return {'denominator': len(values), 'counts': {k:c[k] for k in categories}}


def component_table(pairs):
    local = [classification(x) for x,_ in pairs]
    points = [classification(c.get('point') if c else None) for _,c in pairs]
    bound = [certificate(c.get('bounds') if c else None) for _,c in pairs]
    bp = counts([bypass(x,c) for x,c in pairs], BYPASS)
    nlocal = local.count('joint_only_obstruction')
    bp.update(conditional_local_joint_only_denominator=nlocal,
        conditional_fractions={k:Q(bp['counts'][k],nlocal) if nlocal else None for k in BYPASS[:3]})
    return {'eligible_diagonals':len(pairs), 'local_classifications':counts(local,CLASS),
        'component_point_classifications':counts(points,CLASS),
        'available_component_points':len(points)-points.count('point_unavailable'),
        'unavailable_component_points':points.count('point_unavailable'),
        'complete_component_fields':sum(c is not None and c.get('complete') is True for _,c in pairs),
        'component_bound_statuses':counts(bound,CERT), 'bypass_outcomes':bp}


def write(path, value):
    def convert(item):
        if isinstance(item,Q): return float(item)
        raise TypeError(type(item).__name__)
    with path.open('x',encoding='utf-8') as f:
        json.dump(value,f,indent=2,sort_keys=True,allow_nan=False,default=convert);f.write('\n')


def main():
    inputs = {}
    for role,(rel,pin) in PINS.items():
        path = REPO/rel
        require(sha(path)==pin,'Input hash drift: '+role)
        inputs[role]={'path':rel,'sha256':pin}
    e=json.loads((REPO/PINS['evaluation'][0]).read_bytes())
    d=json.loads((REPO/PINS['diagnostics'][0]).read_bytes())
    s=json.loads((REPO/PINS['split'][0]).read_bytes())
    require(e['arm']==d['arm']=='development','Development arm only')
    require(e['status']=='arm_evaluation_complete' and e['interpretation_permitted'] is True,'Evaluation provenance hold')
    require(d['status']=='supplemental_diagnostics_complete' and d['interpretation_permitted'] is True,'Supplement provenance hold')
    equal(d['input_hashes'],{'evaluation':PINS['evaluation'][1],'split':PINS['split'][1],'diagnostic_script':SCRIPT_SHA},'provenance')
    require(e['source_hashes']['protocol']==PINS['protocol'][1],'Protocol hash agreement')
    require(e['source_hashes']['split']==PINS['split'][1],'Split hash agreement')
    units={u['square_index']:u for u in s['units'] if u['arm']=='development'}
    gf={seed:','.join(g) for g in s['seed_groups'] for seed in g}
    groups=sorted({gf[u['seed']] for u in units.values()},key=lambda x:tuple(map(int,x.split(','))))
    require(len(units)==23 and len(groups)==16,'Frozen selected design size')
    descriptions=e['descriptions']; records=e['prediction_evaluation']['records']
    for name,rows in [('descriptions',descriptions),('prediction',records),('coverage',e['coverage'])]:
        require(len(rows)==23 and {r['square_index'] for r in rows}==set(units),name+' complete inventory')
        for r in rows:
            require(r['seed']==units[r['square_index']]['seed'] and r['group']==gf[r['seed']],name+' identity')
    coverage={r['square_index']:r for r in e['coverage']}
    dm={r['square_index']:r for r in descriptions}
    for row in descriptions:
        u=units[row['square_index']]
        require(row['eligible']==coverage[row['square_index']]['eligible'],'Coverage agreement')
        expected={tuple(x) for x in u['diagonals']} if row['eligible'] else set()
        require(len(row['diagonals'])==(2 if row['eligible'] else 0) and {tuple(x['endpoint_IDs']) for x in row['diagonals']}==expected,'Diagonal inventory')
    scored=[]; labels=[]
    for r in records:
        if r['status']!='scored':continue
        u=units[r['square_index']]; target=u['withheld_ID']; coords=u['coordinate_order_00_01_10_11']
        k=coords.index(target); endpoints=[coords[k^1],coords[k^2]]
        diag=next(z for z in dm[r['square_index']]['diagonals'] if set(z['endpoint_IDs'])==set(endpoints))
        require(r['true_joint_only_obstruction']==(diag['primary']['classification']=='joint_only_obstruction'),'Classification endpoint labels')
        for method in METHODS:
            m=r['methods'][method]
            equal(m['classification_error'],int(m['predicted_joint_only_obstruction']!=r['true_joint_only_obstruction']),'Stored classification loss')
            for j in range(2):
                err=Q(m['predicted_target'][j])-Q(r['observed_target'][j])
                for key,value in [('signed_error',err),('absolute_error',abs(err)),('squared_error',err*err)]:
                    equal(m[key][j],value,'Stored '+key)
        scored.append(dict(r,withheld_ID=target))
        labels.append({'square_index':r['square_index'],'group':r['group'],'withheld_ID':target,'endpoint_IDs':endpoints,
            'true_label':r['true_joint_only_obstruction'],'method_labels':{m:r['methods'][m]['predicted_joint_only_obstruction'] for m in METHODS}})
    equal(d['prediction_case_labels_preserved'],labels,'Preserved case labels')
    equal(d['all_prediction_record_statuses'],[{'square_index':r['square_index'],'status':r['status'],'reason':r.get('reason',r.get('missing_columns'))} for r in sorted(records,key=lambda r:r['square_index'])],'All prediction statuses')
    def summarize(kept):
        dr=[r for r in descriptions if r['group'] in kept]
        return {'descriptive':descriptive(dr,kept),'case_balanced_prediction':prediction(scored,kept,False),
            'target_balanced_prediction':prediction(scored,kept,True)}
    full=summarize(groups)
    equal(d['full_arm'],full,'Full arm independent summary')
    loo=[dict(omitted_group=g,**summarize([x for x in groups if x!=g])) for g in groups]
    equal(d['leave_one_group_out'],loo,'All16 group deletions')
    for item in [full]+loo:
        case={k:v for k,v in item['case_balanced_prediction'].items() if k!='weighting'}
        target={k:v for k,v in item['target_balanced_prediction'].items() if k!='weighting'}
        require(case==target,'Exact rational case/target identity')
    flat=[flattened(item) for item in loo]
    require(all(set(x)==set(flat[0]) for x in flat),'Deletion metric coverage')
    ranges={}
    for key in flat[0]:
        values=[r[key] for r in flat]; known=[x for x in values if x is not None]
        ranges[key]={'minimum':min(known) if known else None,'maximum':max(known) if known else None,
            'defined_rows':len(known),'undefined_rows':16-len(known),'all_rows':16}
    equal(d['all_leave_one_out_minmax_ranges'],ranges,'All deletion extrema and null counts')
    # Use evaluation records directly, never the diagnostic component summaries.
    diagonals=[(r,z) for r in descriptions for z in r['diagonals']]
    cd=d['observed_component_diagnostics']; require(cd['eligible_diagonals']==len(diagonals),'Component diagonal count')
    component_rows=[]; flips=Counter(dict.fromkeys(('local_point','component_point_or_bound','local_range','component_range','scenario_local','scenario_component_point_or_bound'),0))
    def state(c):return (classification(c['point']),certificate(c['bounds']))
    for r,z in diagonals:
        rr={k:z[k] for k in ('local_alternative_tolerance','point_tolerance','observed_component','observed_component_alternative_tolerance','local_range_certificate','component_range_certificate','interval_tolerance','local_range_alternative_tolerance','component_range_alternative_tolerance','native_shortfalls')}
        rr.update({k:r[k] for k in ('square_index','seed','group')})
        rr.update(endpoint_IDs=z['endpoint_IDs'],local_point=z['primary'],
            bypass_status_tolerance_zero=bypass(z['primary'],z['observed_component']),
            bypass_status_alternative_tolerance=bypass(z['local_alternative_tolerance'],z['observed_component_alternative_tolerance']))
        scenarios=z['preparation_scenarios']; require(len(scenarios)==12,'12 scenario inventory')
        require(z['preparation_consistent_joint_obstruction']==all(q['local']['classification']=='joint_only_obstruction' for q in scenarios),'Scenario consistency')
        rr['preparation_scenarios']=[dict(q,bypass_status_tolerance_zero=bypass(q['local'],q['component']),bypass_status_alternative_tolerance=bypass(q['local_alternative_tolerance'],q['component_alternative_tolerance'])) for q in scenarios]
        component_rows.append(rr)
        flips['local_point']+=classification(z['primary'])!=classification(z['local_alternative_tolerance'])
        flips['component_point_or_bound']+=state(z['observed_component'])!=state(z['observed_component_alternative_tolerance'])
        flips['local_range']+=certificate(z['local_range_certificate'])!=certificate(z['local_range_alternative_tolerance'])
        flips['component_range']+=certificate(z['component_range_certificate'])!=certificate(z['component_range_alternative_tolerance'])
        for q in scenarios:
            flips['scenario_local']+=classification(q['local'])!=classification(q['local_alternative_tolerance'])
            flips['scenario_component_point_or_bound']+=state(q['component'])!=state(q['component_alternative_tolerance'])
    equal(cd['records'],component_rows,'Preserved component records')
    equal(cd['flip_counts'],dict(flips),'Six tolerance flip counts')
    ct=cd['count_tables']; require(ct['pooled_across_scenarios'] is False,'No scenario pooling')
    for mode,lk,ck in [('zero_tolerance','primary','observed_component'),('alternative_tolerance','local_alternative_tolerance','observed_component_alternative_tolerance')]:
        equal(ct['primary_point_field'][mode],component_table([(z[lk],z[ck]) for _,z in diagonals]),'Primary component counts/'+mode)
    scenarios_expected={(p,g) for p in ('clean.lin.prot.repOne','clean.lin.prot.repTwo','clean.lin.prot.repThree','clean.lin.prot.repFour') for g in ('growth.R60one','growth.R60two','growth.R60three')}
    require(set(ct['named_preparation_scenarios'])=={p+'|'+g for p,g in scenarios_expected},'Named scenario inventory')
    for p,g in sorted(scenarios_expected):
        scenario_rows=[]
        for _,z in diagonals:
            matches=[q for q in z['preparation_scenarios'] if (q['protein_column'],q['growth_column'])==(p,g)]
            require(len(matches)==1,'Unique named scenario per diagonal')
            scenario_rows.append(matches[0])
        dest=ct['named_preparation_scenarios'][p+'|'+g]
        require(dest['protein_column']==p and dest['growth_column']==g,'Scenario field identity')
        for mode,lk,ck in [('zero_tolerance','local','component'),('alternative_tolerance','local_alternative_tolerance','component_alternative_tolerance')]:
            equal(dest[mode],component_table([(q[lk],q[ck]) for q in scenario_rows]),'Scenario counts/'+p+'/'+g+'/'+mode)
    for mode,lk,ck in [('zero_tolerance','local_range_certificate','component_range_certificate'),('alternative_tolerance','local_range_alternative_tolerance','component_range_alternative_tolerance')]:
        expected={'eligible_diagonals':len(diagonals),'local_certificate_statuses':counts([certificate(z[lk]) for _,z in diagonals],CERT),
            'component_certificate_statuses':counts([certificate(z[ck]) for _,z in diagonals],CERT)}
        equal(ct['empirical_range_field'][mode],expected,'Range counts/'+mode)
    require(all(not x['case_balanced_prediction']['discrimination_screen_passed'] for x in loo),'Actual16 prediction failures')
    interesting=[]
    for r,z in diagonals:
        if z['primary']['classification']=='joint_only_obstruction':
            interesting.append({'square_index':r['square_index'],'group':r['group'],'endpoint_IDs':z['endpoint_IDs'],
                'all12_series_consistent':z['preparation_consistent_joint_obstruction'],
                'range_status':z['local_range_certificate']['status'],
                'observed_component_bypass_status':bypass(z['primary'],z['observed_component']),
                'prediction_square_status':next(x['status'] for x in records if x['square_index']==r['square_index'])})
    for role,(rel,pin) in PINS.items():require(sha(REPO/rel)==pin,'Input changed during review: '+role)
    result={'verdict':'PASS_DEVELOPMENT_SUPPLEMENT_ARITHMETIC','inputs':inputs,'diagnostic_source_sha256':SCRIPT_SHA,
        'reviewer_source_sha256':sha(Path(__file__)),'runtime':platform.python_version(),'checks':dict(CHECKS),'maximum_absolute_numeric_difference':MAX_ABS,
        'comparison_tolerance':{'absolute':1e-12,'relative':1e-12,'scope':'Numerical reproduction only'},
        'all16_deletions_verified':True,'deletion_extrema_paths_verified':len(ranges),'component_fields_verified':28,
        'exact_rational_case_target_identity_full_and_all16_deletions':True,
        'full_arm':full,'leave_one_group_out':loo,'all_leave_one_out_minmax_ranges':ranges,
        'local_joint_obstructions':interesting,'prediction_record_status_counts':dict(Counter(r['status'] for r in records)),
        'all_six_tolerance_flip_counts':dict(flips),
        'scope':'Development saved-output arithmetic only. No raw ZIP or reserved outcomes, author module imports, graph recomputation, fitting or selection. Separately authored review with inherited project context, not blinded.'}
    write(HERE/'independent_results.json',result)
    print(json.dumps({k:result[k] for k in ('verdict','checks','maximum_absolute_numeric_difference','deletion_extrema_paths_verified','local_joint_obstructions','prediction_record_status_counts')},indent=2))


if __name__=='__main__':
    main()
