"""Evaluate the attribution amendment on recorded development worlds only."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.programme.acquisition_benchmark import (audit_prefix, nx_graph,
    paired_summary, truth_status)
from src.acquisition import choose_query

DEVELOPMENT = ROOT/'results/programme/acquisition_v1/development'


def run(output):
    if output.exists():
        raise FileExistsError('The additional baseline requires a new output directory')
    start = perf_counter()
    worlds = json.loads((DEVELOPMENT/'worlds.json').read_text(encoding='utf-8'))
    graphs = json.loads((DEVELOPMENT/'graphs.json').read_text(encoding='utf-8'))
    previous = json.loads((DEVELOPMENT/'policy_results.json').read_text(encoding='utf-8'))
    if len(worlds) != 760 or len(previous) != 3800 or any(
            w['seed'] not in range(1000, 1010) for w in worlds):
        raise ValueError('Expected the unchanged original development cohort')
    graph_map = {g['name']: g for g in graphs}
    records = []
    prefixes = 0
    for world in worlds:
        g = graph_map[world['graph']]
        adjacency = g['adjacency']
        a, b = g['endpoints']
        truth = world['truth']
        actual = truth_status(nx_graph(adjacency), truth, a, b)
        if actual != world['truth_status']:
            raise AssertionError('Recorded truth status differs from independent graph check')
        state = [truth[v] if v in (a, b) or v in world['initial_extra'] else -1
                 for v in range(len(truth))]
        queries = []
        policy_start = perf_counter()
        while True:
            status = audit_prefix(adjacency, state, a, b, actual)
            prefixes += 1
            if status != 'unresolved':
                break
            v = choose_query(adjacency, state, a, b, 'lazy_vertex_bisection',
                             len(queries), None)
            if state[v] != -1 or v in (a, b):
                raise AssertionError('Observed vertex selected')
            state[v] = truth[v]
            queries.append(v)
        if not world['clairvoyant_additional'] <= len(queries) <= len(truth)-2-len(world['initial_extra']):
            raise AssertionError('Certificate bounds violated')
        records.append({k: world[k] for k in ('world_id', 'graph', 'family', 'seed',
            'initial_condition', 'clairvoyant_additional', 'q3_minimax_remaining')} | {
            'policy': 'lazy_vertex_bisection', 'queries': len(queries),
            'query_order': queries, 'certificate': status,
            'policy_seconds': perf_counter()-policy_start})
    summary = paired_summary(previous+records, 'development')
    summary.update({'phase': 'development_additional_baseline',
        'worlds': len(worlds), 'additional_policy_runs': len(records),
        'additional_checked_prefixes': prefixes, 'false_certificates': 0,
        'validation_outcomes_examined': False, 'runtime_seconds': perf_counter()-start,
        'original_receipt_commit': '1e1259477f116ebdc1d69a903bb1552c5c2c92b6',
        'original_receipt_unchanged': True,
        'original_data_sha256': {name: sha256((DEVELOPMENT/name).read_bytes()).hexdigest()
            for name in ('worlds.json', 'graphs.json', 'policy_results.json', 'summary.json')}})
    output.mkdir(parents=True, exist_ok=False)
    (output/'policy_results.json').write_text(json.dumps(records, indent=2)+'\n', encoding='utf-8')
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: summary[k] for k in ('additional_policy_runs', 'additional_checked_prefixes',
        'false_certificates', 'overall', 'efficiency_criterion_met')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
