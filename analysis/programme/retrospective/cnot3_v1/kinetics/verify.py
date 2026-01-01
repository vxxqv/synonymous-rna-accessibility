"""Small exact synthetic audit against an independent four-variable LP oracle.

The oracle starts from k0=b+p, k1=b+delta+c, rL*p<=c<=rU*p.
It does not use the scalar inequalities implemented in kinetic_bounds.py.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

sys.dont_write_bytecode = True
from kinetic_bounds import bound_pathway, fold_and_drop, witness


def solve_square(rows, rhs):
    a = [[Q(x) for x in row] + [Q(value)] for row, value in zip(rows, rhs)]
    for column in range(4):
        pivot = next((i for i in range(column, 4) if a[i][column]), None)
        if pivot is None:
            return None
        a[column], a[pivot] = a[pivot], a[column]
        divisor = a[column][column]
        a[column] = [v / divisor for v in a[column]]
        for i in range(4):
            if i == column:
                continue
            multiplier = a[i][column]
            a[i] = [x - multiplier * y for x, y in zip(a[i], a[column])]
    return tuple(row[4] for row in a)


def polytope_oracle(k0, k1, r, delta, physical):
    l0, u0 = map(Q, k0)
    l1, u1 = map(Q, k1)
    rl, ru = map(Q, r)
    dl, du = map(Q, delta)
    # Variables are p,b,c,delta; every inequality is row*x <= rhs.
    inequalities = [((-1, 0, 0, 0), 0), ((0, -1, 0, 0), 0),
        ((-1, -1, 0, 0), -l0), ((1, 1, 0, 0), u0),
        ((0, -1, -1, -1), -l1), ((0, 1, 1, 1), u1),
        ((rl, 0, -1, 0), 0), ((-ru, 0, 1, 0), 0),
        ((0, 0, 0, -1), -dl), ((0, 0, 0, 1), du)]
    if physical:
        inequalities.append(((0, -1, 0, -1), 0))
    points, systems = set(), 0
    for active in itertools.combinations(inequalities, 4):
        systems += 1
        point = solve_square([row for row, _ in active], [value for _, value in active])
        if point is not None and all(sum(a * b for a, b in zip(row, point)) <= value for row, value in inequalities):
            points.add(point)
    return ((min(p[0] for p in points), max(p[0] for p in points)) if points else None), systems


# Expected intervals below are hand-derived exact answers, not oracle outputs.
# Each fixture gives expected physical and minimal-model projections separately.
FIXTURES = [
    ("complete_removal", (2,2),(1,1),(0,0),(0,0), (1,1),(1,1)),
    ("half_residual", (2,2),(1,1),(Q(1,2),Q(1,2)),(0,0), (2,2),(2,2)),
    ("bounded_residual", (2,2),(1,1),(0,Q(1,2)),(0,0), (1,2),(1,2)),
    ("residual_interval_includes_one", (2,2),(1,1),(0,1),(0,0), (1,2),(1,2)),
    ("no_perturbation", (2,2),(2,2),(1,1),(0,0), (0,2),(0,2)),
    ("unit_residual_contradiction", (2,2),(1,1),(1,1),(0,0), None,None),
    ("unit_residual_signed_shift", (2,2),(1,1),(1,1),(-1,-1), (0,1),(0,2)),
    ("zero_difference_no_shift", (2,2),(2,2),(Q(1,2),Q(1,2)),(0,0), (0,0),(0,0)),
    ("zero_difference_hidden_contribution", (2,2),(2,2),(Q(1,2),Q(1,2)),(0,1), (0,2),(0,2)),
    ("rate_increase_no_shift_rejected", (2,2),(3,3),(Q(1,2),Q(1,2)),(0,0), None,None),
    ("rate_increase_with_shift", (2,2),(3,3),(Q(1,2),Q(1,2)),(1,2), (0,2),(0,2)),
    ("positive_net_drop_may_be_all_background", (2,2),(1,1),(0,0),(-1,0), (0,1),(0,1)),
    ("known_residual_positive_compensation", (3,3),(2,2),(Q(1,2),Q(1,2)),(0,Q(1,2)), (2,3),(2,3)),
    ("physical_background_constraint_matters", (2,2),(Q(1,4),Q(1,4)),(Q(1,2),Q(1,2)),(-2,0), (0,Q(1,2)),(0,2)),
    ("zero_rate_boundary", (0,0),(0,0),(1,1),(0,0), (0,0),(0,0)),
    ("rate_boxes_full_removal", (2,3),(1,2),(0,0),(0,0), (0,2),(0,2)),
    ("rate_boxes_partial_residual", (2,3),(1,2),(Q(1,2),Q(1,2)),(0,0), (0,3),(0,3)),
    ("rate_boxes_unit_residual_overlap", (2,3),(1,2),(1,1),(0,0), (0,2),(0,2)),
    ("rate_boxes_unit_residual_disjoint", (3,4),(1,2),(1,1),(0,0), None,None),
    ("rate_boxes_physical_vs_minimal", (2,3),(0,Q(1,2)),(Q(1,2),Q(1,2)),(-3,0), (0,1),(0,3)),
    ("positive_shift_rate_boxes", (2,3),(1,2),(0,1),(1,2), (1,3),(1,3)),
    ("shift_exceeds_upper_total_rate", (2,3),(1,2),(0,1),(3,4), None,None),
    ("zero_perturbed_rate_boundary", (2,2),(0,0),(0,0),(0,0), (2,2),(2,2)),
    ("incompatible_negative_shift_unit_residual", (1,1),(1,1),(1,1),(-2,-1), None,None),
]


def observed(parameters):
    b, p, r, delta = map(Q, parameters)
    assert b >= 0 and p >= 0 and 0 <= r <= 1 and b + delta >= 0
    return b + p, b + delta + r * p


def run():
    start = perf_counter()
    rows, witness_count, systems = [], 0, 0
    for label, k0, k1, r, delta, physical_expected, minimal_expected in FIXTURES:
        for physical, expected in ((True, physical_expected), (False, minimal_expected)):
            expected = tuple(map(Q, expected)) if expected is not None else None
            result = bound_pathway(k0, k1, r, delta, physical)
            oracle, tested = polytope_oracle(k0, k1, r, delta, physical)
            systems += tested
            assert result["p_interval"] == expected == oracle, (label, physical, result, expected, oracle)
            witnesses = []
            if expected is not None:
                for p in sorted(set([expected[0], expected[1], sum(expected) / 2])):
                    parameters = witness(p, k0, k1, r, delta, physical)
                    witnesses.append(parameters)
                    witness_count += 1
            rows.append({"case":label, "physical_background_nonnegative":physical,
                         "p_interval":result["p_interval"], "independent_polytope_interval":oracle,
                         "constructive_witnesses":witnesses})
    alternatives = [
        ("unknown_residual_same_rates", (1,1,0,0), (0,2,Q(1,2),0)),
        ("known_residual_unknown_positive_shift", (1,2,Q(1,2),0), (0,3,Q(1,2),Q(1,2))),
        ("zero_net_response_hides_pathway", (3,0,Q(1,2),0), (1,2,Q(1,2),1)),
        ("positive_drop_can_be_background_only", (2,0,0,-1), (1,1,0,0)),
    ]
    nonidentifiability = []
    for label, first, second in alternatives:
        assert observed(first) == observed(second) and Q(first[1]) != Q(second[1])
        nonidentifiability.append({"case":label, "observed_rates":observed(first),
                                  "first_b_p_r_delta":first, "second_b_p_r_delta":second})
    small, large = fold_and_drop(1,Q(1,2)), fold_and_drop(10,5)
    assert small["half_life_fold_t1_over_t0"] == large["half_life_fold_t1_over_t0"] == 2
    assert large["absolute_net_rate_reduction"] == 10 * small["absolute_net_rate_reduction"]
    a, b = fold_and_drop(4,2), fold_and_drop(1,Q(1,4))
    assert a["half_life_fold_t1_over_t0"] < b["half_life_fold_t1_over_t0"]
    assert a["absolute_net_rate_reduction"] > b["absolute_net_rate_reduction"]
    assert observed((Q(1,100),Q(9,100),0,0))[0] < observed((1,0,0,0))[0]
    return {"status":"passed", "coefficient_domain":"Exact Fraction; no fitted or measured values",
            "named_fixtures":len(FIXTURES), "model_assumption_cases":len(rows),
            "independent_active_constraint_systems":systems, "constructive_witness_count":witness_count,
            "bounded_cases":rows, "point_identification_counterexamples":nonidentifiability,
            "fold_change_counterexamples":{"same_fold_different_absolute_drop":[small,large],
                                          "rank_reversal":[a,b]},
            "scope":"24 hand-declared rational fixtures, not exhaustive parameter coverage; universal bounds use the separate algebraic derivation",
            "runtime_seconds":perf_counter()-start}


def main():
    if not __debug__:
        raise RuntimeError("Assertions must be enabled")
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Preserve previous verification")
    args.out.mkdir(parents=True)
    here = Path(__file__).resolve().parent
    before = {name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in ("kinetic_bounds.py","verify.py")}
    result = run()
    after = {name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in before}
    assert before == after
    result.update({"source_sha256_before":before,"source_sha256_after":after,
                   "software":{"Python":platform.python_version()},"randomness_used":False})
    (args.out/"results.json").write_text(json.dumps(result,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','named_fixtures','model_assumption_cases','independent_active_constraint_systems','constructive_witness_count','runtime_seconds')}))


if __name__ == '__main__':
    raise SystemExit(main())
