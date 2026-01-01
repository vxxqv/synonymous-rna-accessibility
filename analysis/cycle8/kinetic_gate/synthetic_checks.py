"""Exact invented arithmetic checks. No assay or biological outcome inputs."""
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path


def fiber(n, c, p, e, dn, dp):
    if min(n, c, p, e, dn, dp) <= 0:
        raise ValueError('The model requires strict positivity.')
    return (n*(e+dn), e, dn, e*n/c, dp*p/c, dp)


def equilibrium(parameters):
    s, e, dn, dc, k, dp = parameters
    if min(parameters) <= 0:
        raise ValueError('The model requires strict positivity.')
    n = s/(e+dn)
    c = e*n/dc
    p = k*c/dp
    return n, c, p


def balance(amounts, parameters):
    n, c, p = amounts
    s, e, dn, dc, k, dp = parameters
    return s-(e+dn)*n, e*n-dc*c, k*c-dp*p


def reconstruct(n, c, p, a, dc, dp):
    e = dc*c/n
    dn = a-e
    if dn <= 0:
        raise ValueError('Incompatible with strictly positive nuclear degradation.')
    return n*a, e, dn, dc, dp*p/c, dp


def low_bin_and_mean(points, threshold):
    assert sum(w for w, x in points) == 1
    return sum(w for w, x in points if x < threshold), sum(w*x for w, x in points)


def main():
    lattice = (F(1, 3), F(1), F(3))
    count = 0
    reconstruction_count = 0
    for n, c, p, e, dn, dp in itertools.product(lattice, repeat=6):
        params = fiber(n, c, p, e, dn, dp)
        assert balance((n, c, p), params) == (0, 0, 0)
        assert equilibrium(params) == (n, c, p)
        a = e+dn
        dc = params[3]
        assert reconstruct(n, c, p, a, dc, dp) == params
        assert 0 < dc/a < n/c
        scaled = tuple(7*x/5 for x in params)
        assert equilibrium(scaled) == (n, c, p)
        count += 1
        reconstruction_count += 1
    original = (F(8), F(3), F(1), F(2), F(5), F(3))
    alternatives = [original, (F(8), F(1), F(3), F(2, 3), F(5), F(3)), (F(8), F(3), F(1), F(2), F(10), F(6)), tuple(2*x for x in original)]
    assert all(equilibrium(x) == (2, 3, 5) for x in alternatives)
    assert alternatives[1][1]/alternatives[1][5] != original[1]/original[5]
    assert alternatives[2][5]/sum(alternatives[2][1:3]) != original[5]/sum(original[1:3])
    same_bin = [[(F(1,2), F(1,2)), (F(1,2), F(3,2))], [(F(1,2), F(1,2)), (F(1,2), F(7,2))]]
    assert [low_bin_and_mean(x, F(1)) for x in same_bin] == [(F(1,2), F(1)), (F(1,2), F(2))]
    inverted = [[(F(3,5), F(1,2)), (F(2,5), F(10))], [(F(2,5), F(1,2)), (F(3,5), F(3,2))]]
    assert [low_bin_and_mean(x, F(1)) for x in inverted] == [(F(3,5), F(43,10)), (F(2,5), F(11,10))]
    # Same completion-flux/residence-time product, then larger occupancy with smaller flux.
    assert F(1)*2 == F(2)*1
    assert F(1,2)*6 > F(1)*2 and F(1,2) < 1
    ref = (F(2), F(1), F(1), F(1), F(1), F(1))
    changed = (F(7,2), F(1,6), F(1), F(1), F(1), F(1))
    rn, rc, rp = equilibrium(ref)
    vn, vc, vp = equilibrium(changed)
    assert (rn, rc, rp) == (1, 1, 1) and (vn, vc, vp) == (3, F(1,2), F(1,2))
    assert vn+vc > rn+rc and vp < rp and changed[3:] == ref[3:]
    rejects = 0
    for args in [(F(2), F(3), F(5), F(3), F(2), F(1)), (F(2), F(3), F(5), F(2), F(2), F(1))]:
        # a <= d_C*c/n: equality or negative d_N.
        try:
            reconstruct(*args)
        except ValueError:
            rejects += 1
    assert rejects == 2
    try:
        fiber(F(0), F(1), F(1), F(1), F(1), F(1))
    except ValueError:
        rejects += 1
    assert rejects == 3
    out = Path(__file__).with_name('synthetic_results.json')
    out.write_text(json.dumps({'status':'pass', 'arithmetic':'fractions.Fraction, no tolerances', 'fiber_cases':count, 'inverse_reconstruction_cases':reconstruction_count, 'counterexample_parameter_rows':[[str(x) for x in row] for row in alternatives], 'low_bin_equal_probability_different_means':True, 'low_bin_direction_inversion':True, 'occupancy_flux_confounding':True, 'total_RNA_protein_opposition_without_k_or_dP_change':True, 'boundary_rejections':rejects, 'test_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'scope':'Invented finite regression fixtures supplement the algebraic proof; no empirical data and no claim of independent proof verification.'}, indent=2)+'\n', encoding='utf-8')
    print(out.read_text(encoding='utf-8'))


if __name__ == '__main__':
    if not __debug__:
        raise RuntimeError('Optimized Python is not permitted for assertion-based checks.')
    main()
