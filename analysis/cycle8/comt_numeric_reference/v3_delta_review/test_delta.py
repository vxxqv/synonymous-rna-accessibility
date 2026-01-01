"""Bounded v3 review. Imports candidate normalization only for invented fixtures."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
import argparse
import copy
import csv
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import subprocess
import sys
sys.dont_write_bytecode = True


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def fixtures():
    # (two SNP depths, two MNP per-base depths, two MNP counts)
    settings = [(100,100,90,90,5,5), (80,80,100,100,5,5),
                (100,100,100,100,5,5), (100,100,80,100,4,6),
                (90,90,79,100,5,6), (120,110,100,110,5,5)]
    rng = random.Random(713721)
    groups = {}
    for case, (d1,d2,m1,m2,c1,c2) in enumerate(settings):
        rows = []
        for i,(ref,alt,pnt,dp,cao) in enumerate([('A','C',1,d1,10),('A','G',1,d2,20),
                                              ('AA','CC',1,m1,c1),('AA','CC',2,m2,c2)]):
            rows.append({'Sample':'INVENTED_S','ROI':'ROI_1','Input':'Total_RNA',
                'Fraction':'not_fractionated','Quadrant':'not_sorted','#CHROM':'toy','POS':1,
                'REF':ref,'ALT':alt,'POS_NT':pnt,'REF_NT':'A','ALT_NT':'C' if alt!='G' else 'G',
                'FILTER':'PASS','UP_REF_NT':'A','DOWN_REF_NT':'C','AA_POS':'1',
                'REF_AA':'Q','ALT_AA':'P','LOCATION':'CDS','DP':dp,'CAO':cao,'CAF':cao/dp,
                'record_id':str(case)+':'+str(i)})
        for order in range(4):
            r = copy.deepcopy(rows)
            if order==1:r.reverse()
            elif order>1:rng.shuffle(r)
            groups[f'case{case}_permutation{order}'] = r
    return groups


def main(candidate,pilot,out):
    candidate,pilot,out=Path(candidate).resolve(),Path(pilot).resolve(),Path(out).resolve()
    if out.exists():raise FileExistsError(out)
    manifest=json.loads((candidate/'candidate_manifest_v3.json').read_text())
    for rel,h in manifest['files'].items():
        if sha(candidate/rel)!=h:raise ValueError('Candidate hash mismatch '+rel)
    pinned_before={p.name:sha(p) for p in candidate.iterdir() if p.is_file()}
    norm=load(candidate/'normalize.py','invented_v3_normalize')
    cases=fixtures();out.mkdir(parents=True)
    inp=out/'invented_rows.tsv'
    fields=['fixture']+[k for k in next(iter(cases.values()))[0] if k!='record_id']
    with inp.open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for name,rows in cases.items():
            for row in rows:w.writerow({'fixture':name,**{k:v for k,v in row.items() if k!='record_id'}})
    rscript=Path(json.loads((pilot/'runtime/runtime_manifest.json').read_text())['Rscript'])
    proc=subprocess.run([str(rscript),'--vanilla',str(Path(__file__).parent/'source_fixture.R'),
                         str(pilot),str(inp),str(out/'source_reference.tsv')],capture_output=True,text=True)
    (out/'source_R_execution.log').write_text(proc.stdout+'\n'+proc.stderr,encoding='utf8')
    if proc.returncode:raise RuntimeError('Original-source invented helper failed')
    with (out/'source_reference.tsv').open(newline='',encoding='utf8') as f:
        source=list(csv.DictReader(f,delimiter='\t'))
    expected={}
    checks=0
    for name,rows in cases.items():
        statuses={r['record_id']:{'VAR_ID':':'.join(str(r[k]) for k in ('#CHROM','POS','REF','ALT')),
                  'frequency_pass':True,'globally_rejected':False} for r in rows}
        contract={'members':{'INVENTED_S':{'sample':'INVENTED_S','ROI':1,'trait':'R'}},
                  'targets':[{'ROI':1,'edit_key':'1:A>C','hgvs_nt':'invented'}]}
        got=norm.normalized(rows,statuses,contract)
        line=got['denominator_lineage'][0]
        erows=[r for r in source if r['fixture']==name]
        assert len(erows)==3
        assert [a['VAR_ID'] for a in line['alleles']]==[r['VAR_ID'] for r in erows]
        assert line['selected_depth_variant']==erows[0]['VAR_ID']
        unique=[]
        for r in erows:
            if F(r['DP']) not in unique:unique.append(F(r['DP']))
            a=next(a for a in line['alleles'] if a['VAR_ID']==r['VAR_ID'])
            assert F(a['CAO'])==F(r['CAO']) and F(a['DP'])==F(r['DP'])
            assert F(line['WT'])==F(r['WT_CAO'])
            exact=(F(r['CAO'])+F(1,2))/(F(line['WT'])+F(1,2))
            assert abs(exact-F(r['q']))<=F(1,10**14)
            checks+=4
        assert [F(x) for x in line['ordered_unique_DP']]==unique
        target=got['measurements'][0]
        assert F(target['q'])==(F(10)+F(1,2))/(F(line['WT'])+F(1,2))
        expected[name]=line
        if name.startswith('case1_'):assert F(line['DP'])==80 and F(line['WT'])==45
    for case in range(6):
        first=expected[f'case{case}_permutation0']
        for order in range(1,4):assert expected[f'case{case}_permutation{order}']==first
    rows=copy.deepcopy(cases['case1_permutation0']);rows[0]['CAO']=81
    statuses={r['record_id']:{'VAR_ID':':'.join(str(r[k]) for k in ('#CHROM','POS','REF','ALT')),
                            'frequency_pass':True,'globally_rejected':False} for r in rows}
    try:norm.normalized(rows,statuses,contract)
    except ValueError as e:assert 'Negative reconstructed WT' in str(e)
    else:raise AssertionError('Negative WT permitted')
    after={p.name:sha(p) for p in candidate.iterdir() if p.is_file()}
    assert after==pinned_before
    receipt={'status':'PASS_INVENTED_ORIGINAL_SOURCE_VERSUS_V3_NORMALIZATION',
      'candidate_manifest_sha256':sha(candidate/'candidate_manifest_v3.json'),
      'normalizer_sha256':sha(candidate/'normalize.py'),'fixture_cases':6,'input_permutations_each':4,
      'collapsed_comparisons':checks,'negative_WT_refused':True,'first_lower_than_maximum_control_passed':True,
      'candidate_files_unchanged':True,'real_CAO_CAF_or_downstream_outcomes_accessed':False,
      'code_sha256':sha(__file__),'R_fixture_sha256':sha(Path(__file__).parent/'source_fixture.R')}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--candidate',required=True);p.add_argument('--pilot',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();main(a.candidate,a.pilot,a.out)
