"""Independent completion oracles for reviewed positive-interval implementation."""
if not __debug__:
    raise RuntimeError('Optimized execution is forbidden.')

from fractions import Fraction as F
from pathlib import Path
import hashlib
import itertools
import json
import types

HERE = Path(__file__).resolve().parent
PILOT = HERE.parents[1] / 'cprime_pilot'


def load_core():
    path = PILOT / 'interval_core.py'
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == '9d0bd0f3363e158e4df1cb028c7922b0bf94284bc6169138f894b072871c4cf9'
    mod = types.ModuleType('independent_review_core')
    exec(compile(data, str(path), 'exec'), mod.__dict__)
    return mod


def direct_extrema(pairs):
    candidates = {x for p in pairs for x in p if x is not None and x > 0}
    candidates |= {min(candidates)/2} if candidates else {F(1)}
    domains = [[x for x in sorted(candidates) if x >= l and (u is None or x <= u)] for l, u in pairs]
    realized = [max(y)/min(y) for y in itertools.product(*domains)]
    high = None if len(pairs) > 1 and any(l == 0 or u is None for l, u in pairs) else max(realized)
    return min(realized), high


def allowed_orders(field):
    vertices = sorted(field)
    pairs = [field[v] for v in vertices]
    candidates = {x for p in pairs for x in p if x is not None and x > 0}
    candidates |= {min(candidates)/2} if candidates else {F(1)}
    domains = [[x for x in sorted(candidates) if x >= l and (u is None or x <= u)] for l, u in pairs]
    orders = tuple(itertools.permutations(range(len(vertices))))
    mask = 0
    for values in itertools.product(*domains):
        for i, order in enumerate(orders):
            if all(values[a] <= values[b] for a, b in zip(order, order[1:])):
                mask |= 1 << i
    return mask


def validate_certificate(result, fields):
    vertices = set(next(iter(fields.values())))
    expected = {(v, w, c) for c, d in fields.items() for v in d for w in d if v != w and d[v][1] is not None and d[v][1] < d[w][0]}
    obtained = {(r['lower_state'], r['higher_state'], o['context']) for r in result['edges'] for o in r['origins']}
    assert expected == obtained
    if result['status'] == 'feasible':
        order = result['latent_order']
        assert len(order) == len(vertices) and set(order) == vertices
        for c, d in fields.items():
            y = result['responses'][c]
            for v in vertices:
                l, u = d[v]
                assert 0 < y[v] and l <= y[v] and (u is None or y[v] <= u)
            assert all(y[a] <= y[b] for a, b in zip(order, order[1:]))
    else:
        cycle = result['cycle']
        assert len(cycle) >= 3 and cycle[0] == cycle[-1]
        assert [(e['lower_state'], e['higher_state']) for e in result['cycle_edges']] == list(zip(cycle, cycle[1:]))
        for edge in result['cycle_edges']:
            a, b = edge['lower_state'], edge['higher_state']
            assert edge['origins']
            for origin in edge['origins']:
                c = origin['context']
                assert origin['upper'] == fields[c][a][1] < fields[c][b][0] == origin['lower']


def main():
    core = load_core()
    base = (F(1,2), F(1), F(2), F(3))
    intervals = tuple((l,u) for l in base for u in base if l<=u) + ((F(0),None),(F(0),F(2)),(F(1),None))
    count = 0
    for n in (1,2,3):
        for pairs in itertools.product(intervals, repeat=n):
            expected = direct_extrema(pairs)
            result = core.marker_bounds(pairs)
            assert (result['minimum'],result['maximum']) == expected
            assert result['maximum_attained'] == (expected[1] is not None)
            previous = result
            for gamma in (F(5,4),F(3,2),F(2)):
                current = core.marker_bounds([core.widen(p,gamma) for p in pairs])
                assert current['minimum'] <= previous['minimum']
                assert previous['maximum'] is None or current['maximum'] is None or current['maximum'] >= previous['maximum']
                previous = current
            count += 1
    certificates = 0
    for n in (2,3):
        patterns = list(itertools.product(intervals,repeat=n))
        for i in range(1024):
            a,b = patterns[(i*131+7)%len(patterns)],patterns[(i*239+17)%len(patterns)]
            fields = {'A':dict(zip(range(n),a)),'B':dict(zip(range(n),b))}
            result = core.common_order(fields)
            assert (result['status']=='feasible') == bool(allowed_orders(fields['A']) & allowed_orders(fields['B']))
            validate_certificate(result,fields)
            certificates += 1
    assert core.marker_bounds([(1,10),(4,6)])['maximum'] == 6
    assert core.marker_bounds([(0,None)]) == {'minimum':F(1),'maximum':F(1),'maximum_attained':True}
    rejects = 0
    for pairs in [[],[(0,0)],[(-1,2)],[(2,1)]]:
        try: core.marker_bounds(pairs)
        except ValueError: rejects+=1
    assert rejects == 4
    result = {'status':'pass','core_sha256':hashlib.sha256((PILOT/'interval_core.py').read_bytes()).hexdigest(),
              'exact_completion_patterns':count,'interval_inventory_size':len(intervals),
              'common_order_finite_completion_oracles':certificates,'all_witnesses_and_edge_origins_checked':True,
              'known_zero_tie_missing_distinct_marker_cases':True,'invalid_interval_rejections':rejects,
              'scope':'Invented fixtures only; independent completion enumeration, no production outcome inputs.'}
    (HERE/'math_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
