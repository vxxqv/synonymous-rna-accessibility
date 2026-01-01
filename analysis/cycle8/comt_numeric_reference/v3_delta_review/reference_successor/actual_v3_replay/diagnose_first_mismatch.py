"""Authorized first-mismatch arithmetic only. Does not refit CV or rerun normalization."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def median_witness(rows, convert):
    bysite = Counter(r['site'] for r in rows)
    values = sorted((convert(r['contrasts']['T']['center']), F(1, bysite[r['site']]), r['hgvs_nt']) for r in rows)
    half = sum(w for v, w, k in values) / 2
    cumulative = F()
    for i, (v, w, key) in enumerate(values):
        cumulative += w
        if cumulative >= half:
            selected = [key] if cumulative > half else [key, values[i+1][2]]
            return (v if cumulative > half else (v + values[i+1][0]) / 2), selected, str(half), str(cumulative)
    raise AssertionError('No median')


def main(root, out):
    root, out = Path(root), Path(out)
    if out.exists():
        raise FileExistsError(out)
    paths = {'contrasts': root/'work/research_cycle8/comt_results/actual_v3/normalize/contrasts.json',
        'predictions': root/'work/research_cycle8/comt_results/actual_v3/predict/predictions.json',
        'kernel': root/'work/research_cycle8/comt_pilot/proxy_models.py',
        'prediction_source': root/'work/research_cycle8/comt_pilot/implementation_v3/prediction.py',
        'reference': root/'work/research_cycle8/comt_numeric_reference/v3_delta_review/reference_successor/reference_source_first.py'}
    hashes = {'contrasts':'5654df32df45840415e483a53127a12744887ac9c6b49a18f0a9387bfa04f083',
        'predictions':'a28d2075d078224afd51328056604a145da1f59928486dac8df857c7c2f76bfe',
        'kernel':'8811af9750e813c798122c7e936b6279a0f3b0e53fea1924cb2d633162a34093',
        'prediction_source':'92d0d46fda3b7b3c49c5a4d8780c8b57bcf9c15d9c8095aeca27555f38e2dede',
        'reference':'222bffadbc94710810f9b9ea9c12b594281d4cb7cc83a719abfc6d6fb16143ec'}
    for role, path in paths.items():
        if sha(path) != hashes[role]:
            raise ValueError('Authorized diagnosis input changed: ' + role)
    data = json.loads(paths['contrasts'].read_text(encoding='utf8'))['code']
    target = next(r for r in data if r['hgvs_nt'] == 'c.18C>A')
    fold = target['assessment_fold']
    training = [r for r in data if r['assessment_fold'] != fold and all(r['contrasts'][k]['complete'] for k in ('R','T','L'))]
    decimal, d_ids, half, cum = median_witness(training, lambda x:F(str(x)))
    binary, b_ids, half_b, cum_b = median_witness(training, F)
    # One baseline for one already-authorized fold, independently before retrieving its answer.
    pp = json.loads(paths['predictions'].read_text(encoding='utf8'))['code']['heavy']
    stored = next(p for p in pp['predictions'] if p['hgvs_nt'] == 'c.18C>A')
    actual = F(stored['predictions']['pooled_median']['exact_stored_input_arithmetic'])
    assert stored['fold'] == fold and binary == actual and decimal != actual
    assert d_ids == b_ids and half == half_b and cum == cum_b
    tiny_decimal, tiny_binary = F(str(0.1)), F(0.1)
    assert tiny_decimal == F(1,10) and tiny_binary != tiny_decimal
    result = {'status':'FIRST_MISMATCH_EXPLAINED_BY_REFERENCE_INPUT_RATIONALIZATION',
        'scope':'One pooled-median fit at the already-failed target/fold only; no normalization/CV rerun, no model or output repair',
        'variant':'c.18C>A','map':'code','mode':'heavy','fold':fold,
        'training_rows':len(training),'training_sites':len({r['site'] for r in training}),
        'median_selected_IDs':d_ids,'half_weight':half,'cumulative_weight_at_selection':cum,
        'reference_decimal_string_fraction':str(decimal),'independently_derived_binary_float_fraction':str(binary),
        'production_stored_fraction':str(actual),'production_stored_float':stored['predictions']['pooled_median']['float'],
        'signed_decimal_reference_minus_production':str(decimal-actual),
        'absolute_rational_difference':str(abs(decimal-actual)),
        'absolute_difference_as_float':float(abs(decimal-actual)),
        'difference_after_float_cast':abs(float(decimal)-float(actual)),
        'relative_scaled_rational_difference_as_float':float(abs(decimal-actual)/max(F(1),abs(decimal),abs(actual))),
        'source_semantics':{'production_kernel_number':'Fraction(x), preserving exact stored Python binary-float value',
                           'independent_reference_rational':'Fraction(str(x)), preserving shortest decimal text as a different rational',
                           'production_serialization':'str(Fraction) plus float(Fraction); no intermediate float conversion of the exact fraction'},
        'minimal_synthetic_counterexample':{'input':0.1,'decimal_reference':str(tiny_decimal),
            'production_binary':str(tiny_binary),'binary_minus_decimal':str(tiny_binary-tiny_decimal),
            'absolute_difference':float(abs(tiny_binary-tiny_decimal)),'float_displays_equal':float(tiny_decimal)==float(tiny_binary)},
        'input_hashes':hashes,'diagnostic_code_sha256':sha(__file__),
        'full_prediction_or_metric_agreement_established':False,'frozen_code_or_tolerances_modified':False}
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main(sys.argv[1],sys.argv[2])
