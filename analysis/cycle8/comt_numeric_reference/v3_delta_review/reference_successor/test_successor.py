"""No real outcomes. Source-rule synthetic checks and exact predecessor-function comparison."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
import argparse
import ast
import copy
import csv
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True
import reference_source_first as new


def load(path):
    spec = importlib.util.spec_from_file_location('old_reference_only', path)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    return old


def functions(path):
    text = Path(path).read_text(encoding='utf8')
    return {n.name: ast.get_source_segment(text, n) for n in ast.parse(text).body if isinstance(n, ast.FunctionDef)}


def main(predecessor, source_fixture, invented_stage, out):
    predecessor, source_fixture, invented_stage, out = map(Path, (predecessor, source_fixture, invented_stage, out))
    if out.exists():
        raise FileExistsError(out)
    assert new.sha(predecessor) == new.PREDECESSOR_REFERENCE_SHA256
    aa, bb = functions(predecessor), functions(new.__file__)
    changed = sorted(k for k in aa if aa[k] != bb[k])
    assert changed == ['compare_normalization', 'normalize', 'read_tsv', 'run']
    old = load(predecessor)
    cc = json.loads((invented_stage / 'normalization/contrasts.json').read_text())
    independently_derived = new.downstream(cc)
    assert independently_derived == old.downstream(cc)
    pp = json.loads((invented_stage / 'prediction/predictions.json').read_text())
    ee = json.loads((invented_stage / 'evaluation/evaluation.json').read_text())
    comparison = new.compare_downstream(independently_derived, pp, ee)
    # The preserved source-helper receipt uses invented values only, independently checked earlier.
    with (source_fixture / 'invented_rows.tsv').open(newline='', encoding='utf8') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    with (source_fixture / 'source_reference.tsv').open(newline='', encoding='utf8') as f:
        source = list(csv.DictReader(f, delimiter='\t'))
    checks = 0
    for fixture in sorted({r['fixture'] for r in rows}):
        raw = [{**r, 'record_id': fixture + ':' + str(i)} for i, r in enumerate(rows) if r['fixture'] == fixture]
        masks = [{'record_id': r['record_id'], 'VAR_ID': new.raw_variant_id(r),
                  'frequency_pass': True, 'globally_rejected': False} for r in raw]
        contract = {'members': {'INVENTED_S': {'sample': 'INVENTED_S', 'ROI': 1, 'trait': 'R'}},
                    'targets': [{'ROI': 1, 'edit_key': '1:A>C', 'hgvs_nt': 'invented'}]}
        nn = new.normalize(raw, masks, contract)
        ss = [r for r in source if r['fixture'] == fixture]
        for r in ss:
            c = nn['collapsed'][('INVENTED_S', r['VAR_ID'])]
            assert c['DP'] == F(r['DP']) and c['CAO'] == F(r['CAO']) and c['WT'] == F(r['WT_CAO'])
            assert abs(c['q'] - F(r['q'])) <= F(1, 10**14)
            checks += 4
        if fixture.startswith('case1_'):
            assert nn['collapsed'][('INVENTED_S', 'toy:1:A:C')]['WT'] == 45
    # Exercise the complete new-lineage comparator and its falsifiers on invented first<max.
    raw = [{**r, 'record_id': 'negative:' + str(i)} for i, r in enumerate(rows) if r['fixture'] == 'case1_permutation0']
    masks = [{'record_id': r['record_id'], 'VAR_ID': new.raw_variant_id(r), 'frequency_pass': True, 'globally_rejected': False} for r in raw]
    nn = new.normalize(raw, masks, contract)
    alleles = sorted(next(iter(nn['denominators'].values())), key=lambda a: a['VAR_ID'])
    measure = nn['measurements'][('invented', 'INVENTED_S')]
    stored = {'measurements': [{**measure, 'q': str(measure['q']), 'corrected_q': str(measure['corrected_q'])}],
              'denominator_lineage': [{'Sample': 'INVENTED_S', 'POS': 1, 'DP': '80', 'WT': '45', 'sum_retained_CAO': '35',
                'ordered_unique_DP': ['80', '100'], 'selected_depth_variant': 'toy:1:A:C', 'depth_conflict': True,
                'depth_rule': 'source_sorted_VAR_ID_first_match',
                'alleles': [{'VAR_ID': a['VAR_ID'], 'DP': str(a['DP']), 'CAO': str(a['CAO'])} for a in alleles]}]}
    new.compare_normalization(nn, stored)
    falsifiers = [('DP', '100'), ('selected_depth_variant', 'toy:1:AA:CC'),
                  ('ordered_unique_DP', ['100', '80']), ('depth_rule', 'maximum_DP'), ('depth_conflict', False)]
    for key, value in falsifiers:
        altered = copy.deepcopy(stored)
        altered['denominator_lineage'][0][key] = value
        try:
            new.compare_normalization(nn, altered)
        except AssertionError:
            pass
        else:
            raise AssertionError('Incorrect first-depth lineage accepted: ' + key)
    for r in raw:
        if r['REF'] == 'A' and r['ALT'] == 'C':
            r['CAO'] = '81'
    try:
        new.normalize(raw, masks, contract)
    except new.InvalidDenominatorError as e:
        assert e.details['selected_DP'] == '80' and e.details['WT'] == '-26'
    else:
        raise AssertionError('Invalid denominator repaired or accepted')
    # Reader intentionally does not decode unrelated CAF/RF features.
    out.mkdir(parents=True)
    header = ['record_id','Sample','#CHROM','POS','REF','ALT','POS_NT','REF_NT','ALT_NT','DP','CAO','ROI','CAF']
    poison = out / 'invented_poisoned_unused.tsv'
    record = ['x','INVENTED_S','toy','1','A','C','1','A','C','80','10','ROI_1']
    poison.write_bytes(('\t'.join(header)+'\n'+'\t'.join(record)+'\t').encode()+b'\xffSEALED_CAF\n')
    parsed = new.read_tsv(poison)
    assert parsed[0]['CAO'] == '10' and 'CAF' not in parsed[0]
    for flag in ('-O','-OO'):
        r = subprocess.run([sys.executable,flag,str(Path(new.__file__).resolve()),'--help'],capture_output=True)
        assert r.returncode and b'Optimized Python' in r.stderr
    result = {'status':'PASS_POST_DP_EXPOSURE_PRE_CAO_SOURCE_FIRST_REFERENCE_TESTS',
        'predecessor_function_changes':changed,
        'unchanged_common_functions':sorted(k for k in aa if aa[k] == bb[k]),
        'exact_downstream_objects_equal_to_sealed_predecessor':True,
        'invented_downstream_comparison':comparison,'source_R_synthetic_comparisons':checks,
        'lineage_falsifiers_refused':len(falsifiers),'invalid_denominator_no_repair':True,
        'unused_CAF_bytes_not_decoded':True,'optimized_flags_refused':['-O','-OO'],
        'reference_sha256':new.sha(new.__file__),'test_sha256':new.sha(__file__),
        'actual_CAO_or_downstream_outcomes_opened':False}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--predecessor',required=True);p.add_argument('--source-fixture',required=True);p.add_argument('--invented-stage',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();main(a.predecessor,a.source_fixture,a.invented_stage,a.out)
