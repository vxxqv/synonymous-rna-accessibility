"""Separately authored finite quantifier oracle; synthetic inputs only."""
import sys
if not __debug__:
    raise RuntimeError('Optimized Python refused')
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import importlib.util
import json
import math

BASE=Path(__file__).resolve().parents[1]/'scenario_geometry'
PIN='667422c85a13d5a44761b1248fa0bbae9f2469306b80f546ab9c496e3e7128cf'


def check(ok,msg):
    if not ok:raise AssertionError(msg)


def routes(n,edges,s,t):
    edge_set={frozenset(e) for e in edges}
    interior=[i for i in range(n) if i not in (s,t)]
    result=[]
    for k in range(n-1):
        for middle in permutations(interior,k):
            path=(s,)+middle+(t,)
            if all(frozenset(e) in edge_set for e in zip(path,path[1:])):
                result.append(path)
    return result


def direct(fields,paths,s,t,tau):
    E=min(row[v] for row in fields for v in (s,t))
    if tau>E:return 'ineligible_endpoints'
    if any(all(row[v]>=tau for row in fields for v in path) for path in paths):
        return 'fixed_route'
    accepted=[any(all(row[v]>=tau for v in path) for path in paths) for row in fields]
    if all(accepted):return 'adaptive_routes_only'
    if any(accepted):return 'mixed_scenario_connectivity'
    if any(all(any(row[v]>=tau for row in fields) for v in path) for path in paths):
        return 'joint_obstruction_box_connected'
    return 'upper_box_obstruction'


def metrics(fields,paths,s,t):
    caps=tuple(max((min(row[v] for v in path) for path in paths),default=-math.inf) for row in fields)
    return {'E':min(row[v] for row in fields for v in (s,t)),
            'C_fix':max((min(row[v] for row in fields for v in path) for path in paths),default=-math.inf),
            'C_ad':min(caps),'C_max':max(caps),
            'C_upper':max((min(max(row[v] for row in fields) for v in path) for path in paths),default=-math.inf),
            'scenario_capacities':caps}


def main():
    path=BASE/'scenario_geometry.py'
    check(hashlib.sha256(path.read_bytes()).hexdigest()==PIN,'Held implementation hash drift')
    spec=importlib.util.spec_from_file_location('held_branch_d',path)
    impl=importlib.util.module_from_spec(spec);spec.loader.exec_module(impl)
    counts={'complete_graph_field_cases':0,'pair_metric_records':0,'threshold_checks':0,
            'unknown_cases':0,'completion_pair_checks':0,'bounded_helper_checks':0}
    # All64 graphs on4 vertices and all256 ordered binary two-field assignments.
    n=4;possible=list(combinations(range(n),2))
    for mask in range(1<<len(possible)):
        edges=[e for i,e in enumerate(possible) if mask>>i&1]
        paths={(s,t):routes(n,edges,s,t) for s,t in possible}
        for values in product((Q(0),Q(1)),repeat=2*n):
            fields=[values[:n],values[n:]]
            records=impl.analyze(fields,edges)['pairs']
            check(len(records)==6,'Complete pair coverage')
            counts['complete_graph_field_cases']+=1
            for record in records:
                s,t=record['s'],record['t'];ps=paths[s,t]
                ref=metrics(fields,ps,s,t)
                check(all(record[k]==v for k,v in ref.items()),'Direct metric mismatch')
                counts['pair_metric_records']+=1
                for tau in (Q(0),Q(1,2),Q(1),Q(3,2)):
                    label=direct(fields,ps,s,t,tau)
                    check(impl.classify_exact(record,tau)==label,'Direct quantifier classification mismatch')
                    cells=[x for x in record['intervals'] if x['nonempty'] and
                           (tau>=x['lower'] if x['lower_closed'] else tau>x['lower']) and tau<=x['upper']]
                    check(len(cells)==(label!='ineligible_endpoints'),'Critical interval coverage')
                    check(not cells or cells[0]['label']==label,'Critical interval label')
                    for key,expected in [('all_scenario_obstruction',tau>ref['C_max'] and tau<=ref['E']),
                                         ('upper_box_obstruction',tau>ref['C_upper'] and tau<=ref['E'])]:
                        x=record[key]
                        belongs=x['nonempty'] and (tau>=x['lower'] if x['lower_closed'] else tau>x['lower']) and tau<=x['upper']
                        check(belongs==expected,'Union obstruction interval')
                    counts['threshold_checks']+=1
    # Every incomplete two-field ternary(0,1,None) assignment on a3-path.
    # Unknown cells get values0,1,2, including values above the known range.
    edges=[(0,1),(1,2)];pairs=list(combinations(range(3),2))
    ps={pair:routes(3,edges,*pair) for pair in pairs}
    for flat in product((Q(0),Q(1),None),repeat=6):
        missing=[i for i,v in enumerate(flat) if v is None]
        if not missing:continue
        records=impl.analyze([flat[:3],flat[3:]],edges)['pairs']
        check(all(x['status']=='bounded_unknown' and x['intervals'] is None for x in records),'Unknown exact classification leak')
        counts['unknown_cases']+=1
        for values in product((Q(0),Q(1),Q(2)),repeat=len(missing)):
            filled=list(flat)
            for k,v in zip(missing,values):filled[k]=v
            fields=[filled[:3],filled[3:]]
            for record in records:
                s,t=record['s'],record['t'];ref=metrics(fields,ps[s,t],s,t)
                for key,(lo,hi) in record['bounds'].items():
                    check(lo<=ref[key]<=hi,'Metric bound unsound')
                for actual,(lo,hi) in zip(ref['scenario_capacities'],record['scenario_capacity_bounds']):
                    check(lo<=actual<=hi,'Scenario bound unsound')
                counts['completion_pair_checks']+=1
                for tau in (Q(0),Q(1,2),Q(1),Q(2),Q(3)):
                    label=impl.classify_bounded(record,tau)
                    if not label.startswith('unresolved'):
                        check(label==direct(fields,ps[s,t],s,t,tau),'False bounded certificate')
                    counts['bounded_helper_checks']+=1
    # One scenario collapses fixed/adaptive/mixed gaps; connected zero boundary.
    for fields,edges in [([[Q(0),Q(0)]],[(0,1)]),([[Q(1),Q(2)]],[]),([[Q(1)]],[])]:
        r=impl.analyze(fields,edges)
        check(r['pair_count']==len(fields[0])*(len(fields[0])-1)//2,'Unary/one-field pair coverage')
        for rec in r['pairs']:
            check(rec['C_fix']==rec['C_ad']==rec['C_max']==rec['C_upper'],'Single scenario collapse')
    # Full required coverage with nonbiological constant fields; no outcome access.
    result=impl.analyze([[Q(1)]*72 for _ in range(27)],[(k,k+1) for k in range(71)])
    check(result['pair_count']==2556 and sum(len(r['scenario_capacities']) for r in result['pairs'])==69012,'Full coverage authentication')
    receipt={'status':'PASS_SYNTHETIC_INDEPENDENT_CHALLENGE','counts':counts,
             'full_fixture_pairs':2556,'full_fixture_scenario_capacity_entries':69012,
             'implementation_sha256':PIN,'review_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'oracle':'Permutation enumeration of simple routes and literal quantifier evaluation; no production closure copied.',
             'scope':'Synthetic fields only, not observed successor fields, sampling calibration or biological validation.'}
    out=Path(__file__).with_name('receipt.json')
    if out.exists():raise FileExistsError(out)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))


if __name__=='__main__':main()
