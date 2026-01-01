"""Exact-threshold graph certificates over stored engine-disagreement intervals.

The box is a declared sensitivity set, not a statistical confidence region.
No additional folding calculations are performed.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import sys

sys.path.insert(0, str(Path(__file__).parent))
from design_controls import cube, adjacency, bottleneck

ROOT = Path(__file__).resolve().parents[2]


def labels(graph, active):
    seen, result = set(), [-1]*len(graph)
    group = 0
    for u in range(len(graph)):
        if not active[u] or u in seen:
            continue
        seen.add(u)
        stack = [u]
        result[u] = group
        while stack:
            for v in graph[stack.pop()]:
                if active[v] and v not in seen:
                    seen.add(v)
                    result[v] = group
                    stack.append(v)
        group += 1
    return result


def certificate_scan(graph, lower, upper):
    assert len(graph) == len(lower) == len(upper)
    assert all(0 <= lo <= hi <= 1 for lo,hi in zip(lower,upper))
    baseline = labels(graph, [True]*len(graph))
    critical = sorted(set([0.,1.] + lower + upper))
    spans, witness = [], None
    for lo, hi in zip(critical, critical[1:]):
        t = (lo+hi)/2
        possible = labels(graph, [v >= t for v in upper])
        representatives = {}
        candidate = None
        for i, value in enumerate(lower):
            if value < t:
                continue
            base = baseline[i]
            if base in representatives:
                first = representatives[base]
                if possible[first] != possible[i]:
                    candidate = first, i
                    break
            else:
                representatives[base] = i
        if candidate is not None:
            spans.append([lo,hi])
            a,b = candidate
            barrier_lower = min(lower[a],lower[b])-bottleneck(graph,upper,a,b)
            if witness is None or barrier_lower > witness['certified_lifetime_lower']:
                witness = {'a':a,'b':b,'certified_lifetime_lower':barrier_lower,
                           'lower_endpoint_minimum':min(lower[a],lower[b]),
                           'upper_path_bottleneck':bottleneck(graph,upper,a,b)}
    if witness is not None:
        # Recheck the declared witness with decimal rationals and independent BFS.
        a,b = witness['a'],witness['b']
        rational_upper = [Fraction(str(v)) for v in upper]
        bottom = bottleneck(graph,rational_upper,a,b)
        top = min(Fraction(str(lower[a])),Fraction(str(lower[b])))
        t = (bottom+top)/2
        assert top > bottom
        possible = labels(graph,[v >= t for v in rational_upper])
        assert possible[a] != possible[b]
        witness['rational_barrier_interval_open_closed'] = [str(bottom),str(top)]
    return {'any_certified_positive_barrier':bool(spans),
            'threshold_measure_with_some_certified_disconnected_endpoints':sum(hi-lo for lo,hi in spans),
            'critical_intervals':spans,'witness':witness}


def controls():
    graph = adjacency(cube(2))
    positive = certificate_scan(graph,[.9,0.,0.,.9],[1.,.1,.1,1.])
    ambiguous = certificate_scan(graph,[.4,0.,0.,.4],[1.,.6,.6,1.])
    additive = certificate_scan(graph,[0.,.3,.4,.7],[0.,.3,.4,.7])
    assert positive['any_certified_positive_barrier']
    assert abs(positive['witness']['certified_lifetime_lower']-.8) < 1e-12
    assert not ambiguous['any_certified_positive_barrier']
    assert not additive['any_certified_positive_barrier']
    return {'robust_positive':positive,'unresolved_box':ambiguous,'additive_negative':additive}


def run():
    control_results = controls()
    rows = []
    for path in sorted((ROOT/'results/expansion/rnastructure').glob('*_*.json')):
        rs = json.loads(path.read_text())
        vnpath = ROOT/'results/native'/path.name
        vn = json.loads(vnpath.read_text())
        assert rs['genotypes'] == vn['genotypes']
        assert (rs['id'],rs['condition']) == (vn['id'],vn['condition'])
        graph = [[] for _ in vn['genotypes']]
        for i,j in vn['edges']:
            assert sum(a!=b for a,b in zip(vn['genotypes'][i],vn['genotypes'][j])) == 1
            graph[i].append(j)
            graph[j].append(i)
        # Native edges are independently reconstructed to rule out missing edges.
        expected = adjacency([tuple(s) for s in vn['genotypes']])
        assert [sorted(x) for x in graph] == [sorted(x) for x in expected]
        lower = [max(0.,min(v,lo)) for v,(lo,hi) in zip(vn['pair_retention'],rs['retention_rounding_bounds'])]
        upper = [min(1.,max(v,hi)) for v,(lo,hi) in zip(vn['pair_retention'],rs['retention_rounding_bounds'])]
        assert all(lo <= v <= hi for lo,v,hi in zip(lower,rs['pair_retention'],upper))
        scan = certificate_scan(graph,lower,upper)
        if scan['witness']:
            w = scan['witness']
            w['endpoint_genotypes'] = [vn['genotypes'][w['a']],vn['genotypes'][w['b']]]
        rows.append({'id':rs['id'],'condition':rs['condition'],'transcript_id':rs['transcript_id'],
                     'vertices':len(graph),'engine_input_sha256':{
                         str(path.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest(),
                         str(vnpath.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(vnpath.read_bytes()).hexdigest()},
                     'Vienna_longest_barrier':rs['vienna_filtration']['longest_lifetime'],
                     'RNAstructure_longest_barrier':rs['filtration']['longest_lifetime'],
                     'mean_box_width':sum(hi-lo for lo,hi in zip(lower,upper))/len(lower),
                     **scan})
    return {'scope':'Sensitivity on 48 discovery-cohort graphs; two-engine intervals are not calibrated biological error bounds',
            'controls':control_results,'rows':rows,
            'summary':{'graphs':len(rows),'transcripts':len(set(x['transcript_id'] for x in rows)),
                       'graphs_with_Vienna_barrier':sum(x['Vienna_longest_barrier']>1e-12 for x in rows),
                       'graphs_with_RNAstructure_barrier':sum(x['RNAstructure_longest_barrier']>1e-12 for x in rows),
                       'graphs_with_certified_box_robust_barrier':sum(x['any_certified_positive_barrier'] for x in rows),
                       'maximum_certified_witness_lifetime_lower':max((x['witness']['certified_lifetime_lower'] for x in rows if x['witness']),default=0.)},
            'interpretation':'A positive certificate holds for every score field in the declared box. Failure to certify is inconclusive, not proof of barrier absence. All thresholds were scanned. Witnesses are not selected intervention successes.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out',required=True)
    path=Path(parser.parse_args().out)
    result=run()
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result['summary']))
