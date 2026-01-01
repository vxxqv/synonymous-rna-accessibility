"""D2 exact boundary follow-up; existing candidate solver, independent checking."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import argparse
import hashlib
import importlib.util
import json
import platform
import sys

def need(ok, text):
    if not ok:
        raise RuntimeError(text)

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p, j):
    Path(p).write_text(json.dumps(j, indent=2, default=str)+'\n', encoding='utf-8')

def load_solver(path):
    spec = importlib.util.spec_from_file_location('held_d1_solver', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def rational_system(matrix, target):
    """Exact elimination, returning a free-zero solution and pivot columns."""
    n = len(matrix[0])
    a = [list(map(F, row))+[F(t)] for row, t in zip(matrix, target)]
    pivots = []
    for col in range(n):
        row = next((r for r in range(len(pivots), len(a)) if a[r][col]), None)
        if row is None:
            continue
        k = len(pivots)
        a[k], a[row] = a[row], a[k]
        factor = a[k][col]
        a[k] = [z/factor for z in a[k]]
        for r in range(len(a)):
            if r != k and a[r][col]:
                factor = a[r][col]
                a[r] = [x-factor*y for x, y in zip(a[r], a[k])]
        pivots.append(col)
        if len(pivots) == len(a):
            break
    need(not any(not any(r[:n]) and r[n] for r in a), 'Inconsistent exact lift')
    beta = [F(0)]*n
    for k, col in enumerate(pivots):
        beta[col] = a[k][n]
    need([sum(b*x for b, x in zip(beta, row)) for row in matrix] == list(map(F, target)), 'Exact lift residual')
    return beta, pivots

def forced(x, records, intervals):
    pairs, rows = [], []
    for i, (lo, _) in enumerate(intervals):
        for j, (_, up) in enumerate(intervals):
            if records[i]['cell'] == records[j]['cell'] and up is not None and lo > up:
                pairs.append([i, j])
                rows.append([a-b for a, b in zip(x[i], x[j])])
    return pairs, rows

def make_witness(x, records, intervals, beta):
    pairs, rows = forced(x, records, intervals)
    gaps = [sum(a*b for a, b in zip(row, beta)) for row in rows]
    need(not gaps or min(gaps) > 0, 'Nonpositive lift margin')
    if gaps:
        beta = [b/min(gaps) for b in beta]
    latent = [sum(a*b for a, b in zip(row, beta)) for row in x]
    response = [None]*len(x)
    for cell in sorted({r['cell'] for r in records}):
        inds = [i for i, r in enumerate(records) if r['cell'] == cell]
        positive = [a for i in inds for a in intervals[i] if a is not None and a > 0]
        previous = min(positive)/2 if positive else F(1)
        for q in sorted({latent[i] for i in inds}):
            group = [i for i in inds if latent[i] == q]
            previous = max(previous, max(intervals[i][0] for i in group))
            for i in group:
                response[i] = previous
    result = {'status': 'feasible', 'beta': beta, 'latent': latent, 'responses': response,
              'forced_pairs': pairs, 'margins': [sum(a*b for a, b in zip(row, beta)) for row in rows]}
    check(x, records, intervals, result)
    return result

def check(x, records, intervals, result):
    """Independently reconstruct every constraint and verify exact evidence."""
    pairs, rows = forced(x, records, intervals)
    need(result['forced_pairs'] == pairs, 'Constraint identity/coverage mismatch')
    if result['status'] == 'feasible':
        b, z, y = (list(map(F, result[k])) for k in ('beta', 'latent', 'responses'))
        need(len(b) == len(x[0]) and len(z) == len(x) and len(y) == len(x), 'Witness dimensions')
        need(z == [sum(a*c for a, c in zip(row, b)) for row in x], 'Latent equality')
        gaps = [sum(a*c for a, c in zip(row, b)) for row in rows]
        need(gaps == list(map(F, result['margins'])) and all(g >= 1 for g in gaps), 'Exact unit margins')
        for i, (lo, up) in enumerate(intervals):
            need(y[i] > 0 and lo <= y[i] and (up is None or y[i] <= up), 'Response interval')
            for j in range(len(x)):
                if records[i]['cell'] == records[j]['cell'] and z[i] <= z[j]:
                    need(y[i] <= y[j], 'Response monotonicity, including ties')
    elif result['status'] == 'infeasible':
        terms = result['dual']
        need(len({e['constraint_index'] for e in terms}) == len(terms), 'Duplicate dual entry')
        weights = [F(e['weight']) for e in terms]
        need(all(w > 0 for w in weights) and sum(weights) == 1, 'Positive normalized dual')
        selected = []
        for term in terms:
            k = term['constraint_index']
            need(type(k) is int and 0 <= k < len(rows), 'Dual index')
            need(rows[k] == term['difference'], 'Dual row identity')
            selected.append(rows[k])
        need(all(sum(w*row[c] for w, row in zip(weights, selected)) == 0 for c in range(len(x[0]))), 'Dual cancellation')
    else:
        raise RuntimeError('Unresolved solver status')

def synthetic(solver):
    checks = 0
    records = [{'state':str(i), 'cell':'A', 'exon':1} for i in range(3)]
    x = [[1, 0], [1, 1], [1, 2]]
    intervals = [(F(1), F(1)), (F(3), F(3)), (F(2), F(2))]
    result = solver.solve(x, records, intervals)
    check(x, records, intervals, result)
    need(result['status'] == 'infeasible', 'Rank-deficient additive fixture must reject acyclic order')
    checks += 1
    widened = [(lo, F(3, 2)*up) for lo, up in intervals]
    exact = solver.solve(x, records, widened)
    check(x, records, widened, exact)
    need(exact['status'] == 'feasible' and [1,2] not in exact['forced_pairs'], 'Strict boundary tie fixture')
    checks += 1
    before = solver.solve(x, records, [(lo, F(149,100)*up) for lo, up in intervals])
    check(x, records, [(lo, F(149,100)*up) for lo, up in intervals], before)
    need(before['status'] == 'infeasible', 'Preboundary fixture')
    checks += 1
    full = [[1, 0, 0], [1, 1, 0], [1, 0, 1]]
    beta, pivots = rational_system(full, [0, 2, 1])
    need(len(pivots) == 3, 'Full rank lift')
    make_witness(full, records, intervals, beta)
    checks += 1
    unknown = [(F(0),None)]*3
    empty = solver.solve(x, records, unknown)
    check(x, records, unknown, empty)
    need(empty['status'] == 'feasible' and not empty['forced_pairs'], 'All unknown fixture')
    checks += 1
    tied = solver.solve(x, records, [(F(1),F(1))]*3)
    check(x, records, [(F(1),F(1))]*3, tied)
    checks += 1
    for perm in ([2,0,1], [1,2,0]):
        xx, rr, ii = ([a[k] for k in perm] for a in (x, records, intervals))
        test = solver.solve(xx, rr, ii)
        check(xx, rr, ii, test)
        need(test['status'] == 'infeasible', 'Row permutation fixture')
        checks += 1
    swapped = [list(reversed(row)) for row in x]
    test = solver.solve(swapped, records, intervals)
    check(swapped, records, intervals, test)
    need(test['status'] == 'infeasible', 'Column permutation fixture')
    checks += 1
    import copy
    for original, design, inter, key in [(exact,x,widened,'beta'), (result,x,intervals,'dual')]:
        bad = copy.deepcopy(original)
        if key == 'beta':
            bad['beta'][0] = str(F(bad['beta'][0])+1)
        else:
            bad['dual'][0]['weight'] = '2'
        caught = False
        try:
            check(design, records, inter, bad)
        except RuntimeError:
            caught = True
        need(caught, 'Tamper refusal')
        checks += 1
    return {'status':'passed', 'checks':checks, 'fixture_boundary_t':'3/2',
            'scope':'Invented inputs only; no saved B1 intervals or outcome breakpoints decoded.'}

def old_case(old, exon, model):
    matches = [r for r in old['results'] if r['gamma']=='1' and r['exon']==exon and r['model']==model and r['arm']=='both']
    need(len(matches)==1, 'Existing case identity')
    return matches[0]

def main():
    if not __debug__:
        raise RuntimeError('Optimized Python refused')
    p = argparse.ArgumentParser()
    p.add_argument('--workspace',required=True)
    p.add_argument('--manifest',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--synthetic',action='store_true')
    p.add_argument('--seal')
    p.add_argument('--authorization')
    args = p.parse_args()
    out = Path(args.out)
    need(not out.exists(), 'Output overwrite refused')
    root = Path(args.workspace)
    manifest_path = Path(args.manifest)
    manifest = json.loads(manifest_path.read_text())
    inputs = {name:root/name for name in manifest['inputs']}
    need(all(sha(path)==manifest['inputs'][name] for name,path in inputs.items()), 'Input hash drift')
    here = Path(__file__).resolve().parent
    need(sha(here/'protocol.md')==manifest['protocol_sha256'], 'Protocol hash drift')
    import numpy, scipy
    need((platform.python_version(),numpy.__version__,scipy.__version__)==('3.12.14','2.5.4','1.18.1'), 'Runtime mismatch')
    solver_path = root/'work/research_cycle8/context_design_audit/nested_model/runner.py'
    solver = load_solver(solver_path)
    if args.synthetic:
        result = synthetic(solver)
        result['code_sha256']=sha(__file__)
        result['input_manifest_sha256']=sha(manifest_path)
        out.mkdir()
        save(out/'receipt.json',result)
        print(json.dumps(result))
        return
    need(args.seal is not None and args.authorization is not None, 'Execution seal and authorization required')
    seal = json.loads(Path(args.seal).read_text())
    auth = json.loads(Path(args.authorization).read_text())
    need(seal['code_sha256']==sha(__file__) and seal['protocol_sha256']==sha(here/'protocol.md') and seal['input_manifest_sha256']==sha(manifest_path), 'Pre-execution seal mismatch')
    need(auth['seal_sha256']==sha(args.seal) and auth['authorized'] is True, 'Execution not authorized')
    out.mkdir()
    def checkpoint(j):
        save(out/'progress.json',j)
    try:
        b1 = root/'work/research_cycle8/cprime_results/primary_v1/result.json'
        matrix = root/'work/research_cycle8/context_design_audit/inventory_v1/nested_matrices.json'
        mapping = root/'work/research_cycle8/context_design_audit/exact_v1/mapping.json'
        records, base, matrices, _ = solver.load_contract(b1,matrix,mapping)
        boundary = json.loads((root/'work/research_cycle8/cprime_followup/results_v1/result.json').read_text())['exons']
        old = json.loads((root/'work/research_cycle8/context_design_audit/nested_model/actual_v1/result.json').read_text())
        result = {'status':'in_progress','protocol_sha256':sha(here/'protocol.md'),'code_sha256':sha(__file__),
                  'seal_sha256':sha(args.seal),'authorization_sha256':sha(args.authorization),
                  'state_cells':len(records),'states':len({r['state'] for r in records}), 'exons':[], 'm1_at_one':[]}
        checkpoint(result)
        for exon in (17,25,27,39,42):
            indices = [i for i,r in enumerate(records) if r['exon']==exon]
            rr = [records[i] for i in indices]
            xx = [matrices['strict_shared'][i] for i in indices]
            original = [base[i] for i in indices]
            saved = next(b for b in boundary if b['exon']==exon)
            t = F(saved['t_exact'])
            bounds = [(lo, None if up is None else t*up) for lo,up in original]
            unique = sorted({r['state'] for r in rr})
            interval_map = {(r['state'],r['cell']):iv for r,iv in zip(rr,original)}
            order = saved['feasible_order']
            need(sorted(order)==unique, 'Common order identity mismatch')
            rank = {s:k for k,s in enumerate(order)}
            # Validate every saved response and all common-order edges at t.
            for r, (lo,up) in zip(rr,bounds):
                y = F(saved['scaled_response_witnesses'][r['cell']][r['state']])
                need(y>0 and y>=lo and (up is None or y<=up), 'Saved order completion interval')
                for q in rr:
                    if q['cell']==r['cell'] and rank[r['state']]<=rank[q['state']]:
                        need(y<=F(saved['scaled_response_witnesses'][q['cell']][q['state']]), 'Saved order response monotonicity')
            lower_edges=[]
            cycle=saved['lower_certificate']
            if t>1:
                need(cycle and all(cycle[k]['after']==cycle[(k+1)%len(cycle)]['before'] for k in range(len(cycle))), 'Lower cycle not closed')
                for edge in cycle:
                    ratios={cell:interval_map[edge['after'],cell][0]/interval_map[edge['before'],cell][1] for cell in ('K562','KBM7') if interval_map[edge['before'],cell][1] is not None}
                    weight=max(ratios.values())
                    need(weight==F(edge['weight']) and sorted(cell for cell,value in ratios.items() if value==weight)==sorted(edge['contexts']), 'Lower cycle edge origin/weight')
                    lower_edges.append({'before':edge['before'],'after':edge['after'],'weight':weight,'contexts':edge['contexts']})
                need(min(e['weight'] for e in lower_edges)==t, 'Cycle lower-bound threshold mismatch')
            elif t==1:
                need(not cycle, 'Unexpected lower cycle at domain boundary')
            else:
                raise RuntimeError('Boundary outside domain')
            m1=old_case(old,exon,'context_specific_additive')
            need(m1['row_indices']==indices, 'M1 old row indices')
            check([matrices['context_specific_additive'][i] for i in indices],rr,original,m1)
            need(m1['status']=='feasible','M1 at one unexpected rejection')
            result['m1_at_one'].append({'exon':exon,'row_indices':indices,'certificate':m1})
            entry={'exon':exon,'state_count':len(unique),'row_indices':indices,'arbitrary_order_t':t,'lower_cycle':lower_edges,
                   'at_arbitrary_boundary':None,'fallback_evaluations':[]}
            if exon in (25,39,42):
                need(t==1,'Unexpected nonsingular exon boundary')
                certificate=old_case(old,exon,'strict_shared')
                need(certificate['row_indices']==indices,'M0 old row indices')
                check(xx,rr,bounds,certificate)
                need(certificate['status']=='feasible','Old M0 at one rejected')
                entry.update(method='existing_exact_M0_witness_at_one',at_arbitrary_boundary=certificate)
            elif exon==27:
                first=[k for k,r in enumerate(rr) if r['cell']=='K562']
                need([rr[k]['state'] for k in first]==unique,'Single-context state ordering')
                xone=[xx[k] for k in first]
                desired=[rank[s] for s in unique]
                beta,pivots=rational_system(xone,desired)
                need(len(pivots)==len(unique),'Exon27 not full row rank')
                certificate=make_witness(xx,rr,bounds,beta)
                entry.update(method='exact_full_row_rank_lift',rank=len(pivots),pivot_columns=pivots,
                             lifted_common_order=order,at_arbitrary_boundary=certificate)
            else:
                certificate=solver.solve(xx,rr,bounds)
                check(xx,rr,bounds,certificate)
                entry.update(method='existing_solver_new_exact_boundary_test',at_arbitrary_boundary=certificate)
                if certificate['status']=='infeasible':
                    high=F(25,16)
                    points=sorted({t,high}|{lj/ui for (i,(_,ui)) in enumerate(original) if ui is not None for (k,(lj,_)) in enumerate(original) if rr[i]['cell']==rr[k]['cell'] and t<=lj/ui<=high})
                    entry['fallback_breakpoints']=points
                    low_index,high_index=0,len(points)-1
                    def evaluate(value):
                        intervals=[(lo,None if up is None else value*up) for lo,up in original]
                        cert=solver.solve(xx,rr,intervals)
                        check(xx,rr,intervals,cert)
                        entry['fallback_evaluations'].append({'t':value,'certificate':cert})
                        return cert
                    highcert=evaluate(high)
                    need(highcert['status']=='feasible','Known feasible upper bound fails')
                    while high_index-low_index>1:
                        middle=(low_index+high_index)//2
                        cert=evaluate(points[middle])
                        if cert['status']=='feasible':
                            high_index=middle;highcert=cert
                        else:
                            low_index=middle
                    minimum=points[high_index]
                    prior=points[high_index-1]
                    mid=(prior+minimum)/2
                    lowercert=evaluate(mid)
                    need(lowercert['status']=='infeasible','Fallback chamber not infeasible')
                    for term in lowercert['dual']:
                        i,k=lowercert['forced_pairs'][term['constraint_index']]
                        need(original[i][0]/original[k][1]>=minimum,'Dual fails to persist to boundary')
                    entry.update(M0_t=minimum,first_feasible_certificate=highcert,preceding_chamber=[prior,minimum])
            if 'M0_t' not in entry:
                need(entry['at_arbitrary_boundary']['status']=='feasible','No certified threshold')
                entry['M0_t']=t
            entry['same_boundary']=entry['M0_t']==t
            result['exons'].append(entry)
            checkpoint(result)
        need(all(sha(path)==manifest['inputs'][name] for name,path in inputs.items()),'Inputs changed during execution')
        result.update(status='COMPLETE_EXPLORATORY_EXACT_BOUNDARY_AUDIT',all_five_boundaries_equal=all(e['same_boundary'] for e in result['exons']),
                      optimizer_calls=sum(1+len(e['fallback_evaluations']) for e in result['exons'] if e['exon']==17),
                      scope='Exposed descriptive interval compatibility; no biological calibration, new source validation, or generic algorithm priority.')
        save(out/'result.json',result)
        save(out/'receipt.json',{'status':'passed','result_sha256':sha(out/'result.json'),'code_sha256':sha(__file__),
                                 'input_manifest_sha256':sha(manifest_path),'protocol_sha256':sha(here/'protocol.md'),
                                 'seal_sha256':sha(args.seal),'authorization_sha256':sha(args.authorization)})
        print(json.dumps({'status':result['status'],'all_five_boundaries_equal':result['all_five_boundaries_equal'],'optimizer_calls':result['optimizer_calls']}))
    except BaseException as exc:
        save(out/'failure.json',{'status':'STOPPED_UNRESOLVED','type':type(exc).__name__,'message':str(exc)})
        raise

if __name__=='__main__':
    main()
