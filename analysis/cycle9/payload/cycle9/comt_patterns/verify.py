"""Separate direct reconstruction; no production implementation import.

This is a distinct arithmetic formulation by the same project agent. A fresh
reviewer can run/inspect it; it is not itself blinded external verification.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def independent_pair(pair, source):
    maps, allr, allt = {}, [], []
    for mapping in ('code', 'geo'):
        rows = {r['hgvs_nt']:r for r in source[mapping]}
        a, b = rows[pair['a']], rows[pair['b']]
        z = {}; complete = True
        for assay in ('D','R','F'):
            ag = {g['bio']:g for g in a['assays'][assay]['groups']}
            bg = {g['bio']:g for g in b['assays'][assay]['groups']}
            zz = []
            for n in range(1,5):
                x, y = ag.get(n), bg.get(n)
                if x is None or y is None or x['median'] is None or y['median'] is None:
                    zz.append(None); complete = False
                elif x['samples'] != y['samples'] or [v is not None for v in x['values']] != [v is not None for v in y['values']]:
                    zz.append(None); complete = False
                else: zz.append(F(x['median']) - F(y['median']))
            z[assay] = zz
        cases = []
        for i,j,k in itertools.product(range(4),repeat=3):
            R = z['R'][i] - z['D'][j] if z['R'][i] is not None and z['D'][j] is not None else None
            T = z['D'][j] - z['F'][k] if z['D'][j] is not None and z['F'][k] is not None else None
            cases.append((i+1,j+1,k+1,R,T))
        if complete:
            # Independently use separability for the range, not minima of primary records.
            bounds = (min(z['R'])-max(z['D']), max(z['R'])-min(z['D']),
                      min(z['D'])-max(z['F']), max(z['D'])-min(z['F']))
            rsign = {((r>0)-(r<0)) for _,_,_,r,_ in cases}
            tsign = {((t>0)-(t<0)) for _,_,_,_,t in cases}
            category = signs_category(rsign, tsign)
            allr += [r for _,_,_,r,_ in cases]; allt += [t for _,_,_,_,t in cases]
        else: bounds = None; category = 'unresolved_missing_or_misaligned'
        maps[mapping] = {'z':z, 'cases':cases, 'bounds':bounds, 'category':category, 'complete':complete}
    complete = all(x['complete'] for x in maps.values())
    if complete:
        category = signs_category({(x>0)-(x<0) for x in allr}, {(x>0)-(x<0) for x in allt})
        bounds = (min(allr),max(allr),min(allt),max(allt))
        radius = min(abs(x) for x in allr+allt)/4 if category.startswith('reversal_') else F(0)
    else: category = 'unresolved_missing_or_misaligned'; radius = None; bounds = None
    return {'maps':maps,'category':category,'radius':radius,'bounds':bounds,'complete':complete}

def signs_category(R,T):
    fixed = {((1,),(-1,)):'reversal_R_up_T_down', ((-1,),(1,)):'reversal_R_down_T_up',
             ((1,),(1,)):'same_direction_a_higher', ((-1,),(-1,)):'same_direction_b_higher'}
    key = (tuple(sorted(R)),tuple(sorted(T)))
    return fixed.get(key, 'unresolved_with_ties' if 0 in R or 0 in T else 'unresolved_no_ties')

def check(source, result, inventory):
    checks = 0; cases = 0
    def equal(a,b):
        nonlocal checks
        checks += 1
        if a != b: raise ValueError('Independent mismatch: '+repr(a)[:180]+' != '+repr(b)[:180])
    equal([(p['a'],p['b']) for p in result['pairs']],[(p['a'],p['b']) for p in inventory['pairs']])
    names = ('R_lower','R_upper','T_lower','T_upper')
    for row in result['pairs']:
        oracle = independent_pair(row,source)
        equal(row['category'],oracle['category']); equal(row['complete'],oracle['complete'])
        equal(row['unstructured_log_radius_exclusive'],None if oracle['radius'] is None else str(oracle['radius']))
        equal(row['combined_bounds'],None if oracle['bounds'] is None else dict(zip(names,map(str,oracle['bounds']))))
        for mp, ref in oracle['maps'].items():
            actual = row['maps'][mp]
            equal(actual['complete'],ref['complete']);equal(actual['category'],ref['category'])
            equal(actual['same_assay_differences'],{a:[None if x is None else str(x) for x in z] for a,z in ref['z'].items()})
            equal(actual['bounds'],None if ref['bounds'] is None else dict(zip(names,map(str,ref['bounds']))))
            equal(len(actual['joint_cases']),64)
            for rec,(i,j,k,R,T) in zip(actual['joint_cases'],ref['cases']):
                equal(rec,{'R_group':i,'D_group':j,'F_group':k,'DeltaR':None if R is None else str(R),'DeltaT':None if T is None else str(T)})
                cases += 1
            if ref['complete']:
                for name,value in zip(names,ref['bounds']):
                    values=[c[3] if name.startswith('R') else c[4] for c in ref['cases']]
                    equal(actual['extremum_case_indices_zero_based'][name],[j for j,x in enumerate(values) if x==value])
        equal(row['map_categories_differ'], len({r['category'] for r in oracle['maps'].values()}) != 1)
    for name in ('all','SNV'):
        rows = [p for p in result['pairs'] if name=='all' or p['both_SNV']]
        s = result['summaries'][name]
        equal(s['pairs'],len(rows));equal(s['joint_case_records'],128*len(rows))
        equal(s['category_counts'],dict(sorted(Counter(r['category'] for r in rows).items())))
        equal(s['uniform_reversal_pairs'],sum(r['category'].startswith('reversal_') for r in rows))
        fractions=[]
        for site in s['sites_detail']:
            group=[r for r in rows if r['site']==site['site']]
            nr=sum(r['category'].startswith('reversal_') for r in group)
            equal(site['pair_count'],len(group));equal(site['uniform_reversals'],nr)
            equal(site['category_counts'],dict(sorted(Counter(r['category'] for r in group).items())))
            if group:
                rate=F(nr,len(group));fractions.append(rate);equal(site['reversal_fraction'],str(rate))
            else: equal(site['reversal_fraction'],None)
        equal(s['equal_informative_site_reversal_fraction'],str(sum(fractions,F(0))/len(fractions)) if fractions else None)
    return {'checks':checks,'joint_cases_verified':cases,'maximum_exact_discrepancy':'0'}

def main():
    if not __debug__: raise SystemExit('Optimized Python unsupported')
    ap=argparse.ArgumentParser()
    for n in ('input','result','inventory','freeze','out'):ap.add_argument('--'+n,required=True)
    a=ap.parse_args();out=Path(a.out)
    if out.exists():raise FileExistsError(out)
    freeze=json.loads(Path(a.freeze).read_text())
    for n,p in [('input',a.input),('inventory',a.inventory)]:
        if sha(p)!=freeze['input_hashes'][n]:raise ValueError('Frozen input mismatch')
    if sha(__file__)!=freeze['owned_files']['verify.py']:raise ValueError('Verifier source changed')
    paths={n:getattr(a,n) for n in ('input','result','inventory','freeze')}
    before={n:sha(p) for n,p in paths.items()}
    summary=check(json.loads(Path(a.input).read_text()),json.loads(Path(a.result).read_text()),json.loads(Path(a.inventory).read_text()))
    if before!={n:sha(p) for n,p in paths.items()}:raise ValueError('Inputs changed')
    out.mkdir(parents=True)
    payload={'status':'PASS_DISTINCT_DIRECT_RECONSTRUCTION','scope':'Same-author independent arithmetic formulation, not an external biological validation.',
             **summary,'input_hashes':before,'verifier_sha256':sha(__file__)}
    (out/'receipt.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(payload,indent=2))

if __name__=='__main__':main()
