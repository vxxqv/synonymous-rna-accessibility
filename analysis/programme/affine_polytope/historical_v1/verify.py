"""Bounded, reproducible verification of three candidate graph claims.

Production is exact rational. The independent path/cut oracle uses floating
linear programming and an explicit 1e-8 absolute tolerance. It does not use
the production arrangement or widest-path routine. This author check awaits
a fresh independent proof review and reference implementation.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import json
import platform
import random
import time
import argparse
import numpy as np
import scipy
from scipy.optimize import linprog
import networkx as nx
import affine_polytope as ap

HERE = Path(__file__).resolve().parent


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k):encode(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_cases():
    cases = []
    diamond = {'n':4, 'edges':[(0,1),(1,3),(0,2),(2,3)], 's':0, 't':3}
    cases.append({'id':'two_dimensional_corner_failure', **diamond, 'box':[[-1,1],[-1,1]],
                  'coefficients':[[2,0,0],[0,1,1],[0,-1,-1],[2,0,0]]})
    cases.append({'id':'lower_dimensional_line', **diamond, 'box':[[-1,1],[-1,1]],
                  'inequalities':[[[1,-1],0],[[-1,1],0]],
                  'coefficients':[[1,0,0],[0,1,0],[0,0,-1],[1,0,0]]})
    cases.append({'id':'point_parameter_set', **diamond, 'box':[[0,0],[0,0]],
                  'coefficients':[[1,0,0],[0,1,0],[0,0,-1],[1,0,0]]})
    cases.append({'id':'zero_dimension', **diamond, 'box':[],
                  'coefficients':[[1],[0],[-1],[1]]})
    cases.append({'id':'empty_polytope', **diamond, 'box':[[-1,1]],
                  'inequalities':[[[1],-2]], 'coefficients':[[1,0],[0,1],[0,-1],[1,0]]})
    cases.append({'id':'disconnected', 'n':4,'edges':[(0,1),(2,3)],'s':0,'t':3,
                  'box':[[-1,1]],'coefficients':[[1,0],[0,1],[0,-1],[1,0]]})
    cases.append({'id':'same_endpoint', **diamond, 't':0, 'box':[[-1,1],[0,2]],
                  'coefficients':[[1,1,-1],[0,1,0],[0,-1,0],[1,0,0]]})
    cases.append({'id':'direct_edge', **diamond, 'edges':diamond['edges']+[(0,3)],
                  'box':[[-1,1],[-1,1]], 'coefficients':[[1,1,0],[-2,0,1],[0,-1,-1],[1,0,1]]})
    cases.append({'id':'all_tied_negative', **diamond, 'box':[[-1,1],[-1,1]],
                  'coefficients':[[-3,1,-1]]*4})
    rng = random.Random(172904)
    for index in range(72):
        n = 3 + index%3
        d = 1 + (index//3)%3
        edges = {(i,i+1) for i in range(n-1)}
        edges |= {e for e in combinations(range(n),2) if rng.random() < .28}
        coefficients = [[rng.randint(-3,3) for _ in range(d+1)] for _ in range(n)]
        box = [[-1,1] for _ in range(d)]
        inequalities = []
        if index%4 == 0:
            inequalities = [[[1]*d, '1/2']]
        if index%11 == 0:
            box[-1] = [0,0]
        if index%13 == 0:
            coefficients[-1] = coefficients[0].copy()
        cases.append({'id':f'fixed_seed_{index:02d}', 'n':n,'edges':sorted(edges),'s':0,'t':n-1,
                      'box':box,'inequalities':inequalities,'coefficients':coefficients})
    return cases


def path_cut_lp(case):
    """Use Boolean path/cut duality and LPs, not a score arrangement."""
    n, s, t = case['n'], case['s'], case['t']
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    graph.add_edges_from(case['edges'])
    if not nx.has_path(graph,s,t):
        return {'status':'disconnected'}
    coeff = np.array([[float(F(x)) for x in row] for row in case['coefficients']])
    d = len(case['box'])
    if d == 0:
        paths = [[s]] if s == t else list(nx.all_simple_paths(graph,s,t))
        capacity = max(min(coeff[v,0] for v in path) for path in paths)
        barrier = min(coeff[s,0],coeff[t,0])-capacity
        return {'status':'connected','capacity_minimum':capacity,'capacity_maximum':capacity,
                'barrier_minimum':barrier,'barrier_maximum':barrier,'lp_count':0}
    paths = [[s]] if s == t else list(nx.all_simple_paths(graph,s,t))
    cuts = []
    for size in range(1,n+1):
        for selected in combinations(range(n),size):
            group = set(selected)
            if any(set(cut) <= group for cut in cuts):
                continue
            if s in group or t in group or not nx.has_path(graph.subgraph(set(range(n))-group),s,t):
                cuts.append(selected)
    bounds = [(float(F(lo)),float(F(hi))) for lo,hi in case['box']]
    base_A = [[float(F(x)) for x in row] for row,rhs in case.get('inequalities',[])]
    base_b = [float(F(rhs)) for row,rhs in case.get('inequalities',[])]
    lp_count = 0
    def lp(objective, rows, rhs, extra=1):
        nonlocal lp_count
        lp_count += 1
        aa = [row+[0.]*extra for row in base_A] + rows
        bb = base_b + rhs
        result = linprog(objective,A_ub=aa or None,b_ub=bb or None,
                         bounds=bounds+[(None,None)]*extra,method='highs')
        if result.status == 2:
            return None
        if not result.success:
            raise AssertionError((case['id'],result.message))
        return result.fun
    cap_max = []
    bar_min = []
    for path in paths:
        # maximize z with z <= f_v for every path vertex
        value = lp([0.]*d+[-1.], [list(-coeff[v,1:])+[1.] for v in path],
                   [coeff[v,0] for v in path])
        if value is not None:
            cap_max.append(-value)
        for endpoint in {s,t}:
            # minimize z >= f_endpoint - f_v, then minimize over endpoints
            value = lp([0.]*d+[1.], [list(coeff[endpoint,1:]-coeff[v,1:])+[-1.] for v in path],
                       [coeff[v,0]-coeff[endpoint,0] for v in path])
            if value is not None:
                bar_min.append(value)
    cap_min = []
    bar_max = []
    for cut in cuts:
        # minimize z >= f_v for every vertex in the separator
        value = lp([0.]*d+[1.], [list(coeff[v,1:])+[-1.] for v in cut],
                   [-coeff[v,0] for v in cut])
        if value is not None:
            cap_min.append(value)
        # maximize y-z, y <= both endpoint scores and z >= all cut scores
        rows = [list(-coeff[v,1:])+[1.,0.] for v in {s,t}]
        rhs = [coeff[v,0] for v in {s,t}]
        rows += [list(coeff[v,1:])+[0.,-1.] for v in cut]
        rhs += [-coeff[v,0] for v in cut]
        value = lp([0.]*d+[-1.,1.],rows,rhs,extra=2)
        if value is not None:
            bar_max.append(-value)
    if not cap_max:
        return {'status':'empty_parameter_set','lp_count':lp_count}
    return {'status':'connected','capacity_minimum':min(cap_min),'capacity_maximum':max(cap_max),
            'barrier_minimum':min(bar_min),'barrier_maximum':max(bar_max),
            'lp_count':lp_count,'paths':len(paths),'minimal_cuts':len(cuts)}


def classification_check():
    total = connected = unique = multiple = 0
    examples = []
    for n in range(2,6):
        edge_list = list(combinations(range(n),2))
        for bits in range(1 << len(edge_list)):
            total += 1
            graph = nx.Graph()
            graph.add_nodes_from(range(n))
            graph.add_edges_from(e for j,e in enumerate(edge_list) if bits & (1<<j))
            if not nx.has_path(graph,0,n-1):
                continue
            connected += 1
            sets = {frozenset(p) for p in nx.all_simple_paths(graph,0,n-1)}
            minimal = sorted((p for p in sets if not any(q < p for q in sets)), key=lambda p:tuple(sorted(p)))
            if len(minimal) == 1:
                unique += 1
            else:
                multiple += 1
                p,q = minimal[:2]
                common = p & q
                assert not (0 in common and n-1 in common and nx.has_path(graph.subgraph(common),0,n-1))
                # Evaluate the claimed full affine certificate by exact path minima.
                for theta in (F(-1),F(0),F(1)):
                    scores = [F(1) if v in common else theta if v in p-q else -theta if v in q-p else F(-1) for v in range(n)]
                    capacity = max(min(scores[v] for v in path) for path in sets)
                    assert capacity == abs(theta)
                low = [1 if v in common else -1 for v in range(n)]
                assert max(min(low[v] for v in path) for path in sets) == -1
                if len(examples) < 3:
                    examples.append({'n':n,'edges':list(graph.edges),'minimal_sets':[sorted(x) for x in minimal]})
    return {'labelled_graphs':total,'connected_pairs':connected,'unique_minimal_set':unique,
            'multiple_minimal_sets_with_exact_gap':multiple,'examples':examples,
            'scope':'all labelled simple graphs n=2..5 and distinguished endpoints 0,n-1; finite corroboration only'}


def reductions_check():
    clauses = list(product((False,True),repeat=3))
    formulas = satisfiable = unsatisfiable = 0
    witness_checks = 0
    for mask in range(1,1<<8):
        selected = [cl for i,cl in enumerate(clauses) if mask & (1<<i)]
        sat = any(all(any(value == sign for value,sign in zip(assignment,clause)) for clause in selected)
                  for assignment in product((False,True),repeat=3))
        formulas += 1
        satisfiable += int(sat)
        unsatisfiable += int(not sat)
        # Actual undirected graphs for the two constructions.
        n = 2 + 3*len(selected)
        t = n-1
        layered = nx.Graph(); parallel = nx.Graph()
        layered.add_nodes_from(range(n)); parallel.add_nodes_from(range(n))
        groups = [list(range(1+3*i,4+3*i)) for i in range(len(selected))]
        layered.add_edges_from((0,v) for v in groups[0])
        layered.add_edges_from((v,t) for v in groups[-1])
        for before,after in zip(groups,groups[1:]):
            layered.add_edges_from(product(before,after))
        for group in groups:
            parallel.add_edges_from(zip([0]+group,group+[t]))
        possible = False
        guaranteed_grid = True
        # {0,1/2,1} probes ties and fractional witnesses; not proof over continuum.
        for theta in product((F(0),F(1,2),F(1)),repeat=3):
            literal = [F(1)] + [theta[i] if sign else 1-theta[i] for cl in selected for i,sign in enumerate(cl)] + [F(1)]
            complemented = [F(1)] + [1-x for x in literal[1:-1]] + [F(1)]
            possible |= nx.has_path(layered.subgraph([v for v,x in enumerate(literal) if x >= 1]),0,t)
            guaranteed_grid &= nx.has_path(parallel.subgraph([v for v,x in enumerate(complemented) if x >= F(1,2)]),0,t)
            witness_checks += 2
        assert possible == sat
        assert guaranteed_grid == (not sat)
    return {'formulas':formulas,'satisfiable':satisfiable,'unsatisfiable':unsatisfiable,
            'threshold_graph_evaluations':witness_checks,
            'scope':'all nonempty subsets of eight distinct full three-variable clauses; 27 rational parameter points each; proof supplies continuum statement'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if (output/'results.json').exists():
        raise FileExistsError('Refusing to replace a prior result')
    started = time.perf_counter()
    cases = make_cases()
    # Freeze input suite and code before running any new outcome calculation.
    suite = HERE/'cases.json'
    assert json.loads(suite.read_text(encoding='utf-8')) == json.loads(json.dumps(cases)), 'Frozen suite mismatch'
    tracked = [HERE/'claims.md',HERE/'protocol.json',HERE/'affine_polytope.py',Path(__file__),suite]
    inputs = {path.name:sha(path) for path in tracked}
    (output/'freeze.json').write_text(json.dumps({'inputs_sha256':inputs,'suite_count':len(cases),
                    'meaning':'algorithm verification freeze, not independent biological preregistration'},indent=2)+'\n',encoding='utf-8')
    comparisons = []
    for case in cases:
        exact = ap.analyze(case)
        reference = path_cut_lp(case)
        assert exact['status'] == reference['status'], (case['id'],exact,reference)
        errors = {}
        if exact['status'] == 'connected':
            for key in ('capacity_minimum','capacity_maximum','barrier_minimum','barrier_maximum'):
                errors[key] = abs(float(exact[key])-reference[key])
                assert errors[key] <= 1e-8, (case['id'],key,exact,reference)
        comparisons.append({'id':case['id'],'exact':exact,'reference':reference,'absolute_errors':errors})
    classification = classification_check()
    reductions = reductions_check()
    assert inputs == {path.name:sha(path) for path in tracked}, 'Input changed during run'
    result = {'status':'completed','production':'exact rational','lp_reference_tolerance':1e-8,
              'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'networkx':nx.__version__},
              'inputs_sha256':inputs,'cases':comparisons,'classification':classification,'reductions':reductions,
              'seconds':time.perf_counter()-started,'independence':'self authored distinct oracle; fresh independent review pending'}
    (output/'results.json').write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'cases':len(cases),'lp_solves':sum(x['reference'].get('lp_count',0) for x in comparisons),
                      'classification':classification,'reductions':reductions,'seconds':result['seconds']}))


if __name__ == '__main__':
    main()
