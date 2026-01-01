"""One no-Git default-path replay, using authorized saved inputs only."""
import sys
if not __debug__:raise RuntimeError('Optimized checker execution refused')
from pathlib import Path
import gzip,hashlib,json,shutil,subprocess


def main():
    here=Path(__file__).resolve().parent
    workspace=here.parents[2]
    repo=workspace/'outputs/synonymous-rna-accessibility'
    fixture=here/'portable_fixture'
    if fixture.exists():raise FileExistsError(fixture)
    archive=fixture/'archive';out=fixture/'replayed'
    freeze_source=repo/'analysis/cycle7/freeze.json'
    freeze=json.loads(freeze_source.read_text(encoding='utf-8'))
    target=archive/'analysis/cycle7';target.mkdir(parents=True)
    shutil.copyfile(freeze_source,target/'freeze.json')
    for name in freeze['files']:
        dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(freeze_source.parent/name,dest)
    for cycle,name,source in [('cycle6','full_result_v1.json.gz',workspace/'work/research_cycle6/fas_results/full_result_v1.json'),
                              ('cycle7','scenario_geometry.json.gz',workspace/'work/research_cycle7/results/scenario_geometry.json')]:
        dest=archive/'results'/cycle/name;dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_bytes(gzip.compress(source.read_bytes(),mtime=0))
    code=archive/'analysis/cycle7/scenario_independent/verify_observed.py'
    code.parent.mkdir();shutil.copyfile(here/'verify_observed.py',code)
    assert not (archive/'.git').exists()
    cmd=[sys.executable,str(code),'--repo',str(archive),'--out',str(out)]
    p=subprocess.run(cmd,capture_output=True,text=True)
    assert p.returncode==0,p.stderr
    baseline=json.loads((here/'observed_run/verification.json').read_text(encoding='utf-8'))
    actual=json.loads((out/'verification.json').read_text(encoding='utf-8'))
    assert {k:v for k,v in actual.items() if k!='inputs'}=={k:v for k,v in baseline.items() if k!='inputs'}
    assert all(actual['inputs'][k]['decompressed_sha256']==baseline['inputs'][k]['decompressed_sha256'] for k in actual['inputs'])
    assert (out/'all_pairs.csv').read_bytes()==(here/'observed_run/all_pairs.csv').read_bytes()
    assert actual['historical_git_verified'] is False
    checks=[]
    def refusal(arguments,expected):
        p=subprocess.run([sys.executable,str(code),*arguments],capture_output=True,text=True)
        assert p.returncode and expected in p.stderr,(expected,p.stderr)
        checks.append(expected)
    refusal(['--repo',str(archive),'--out',str(out)],'Refusing existing output')
    refusal(['--repo',str(archive),'--out',str(archive/'MUST_NOT_WRITE')],'Output must be outside repository')
    refusal(['--repo',str(fixture/'MISSING_ROOT'),'--out',str(fixture/'MUST_NOT_WRITE')],'Repository directory missing')
    bad=fixture/'altered_freeze.json';bad.write_bytes(b'not JSON, invented altered metadata')
    refusal(['--repo',str(archive),'--out',str(fixture/'MUST_NOT_WRITE'),'--freeze',str(bad)],'Input decompressed hash mismatch')
    assert not (fixture/'MUST_NOT_WRITE').exists() and not (archive/'MUST_NOT_WRITE').exists()
    receipt={'status':'PASS_NO_GIT_PORTABLE_DEFAULT_REPLAY','scientific_verification_identical':True,
             'all_pairs_csv_byte_identical':True,'all_decompressed_input_hashes_identical':True,
             'historical_git_verified':False,'stored_compression_hashes_may_differ':True,
             'guard_checks':checks,'copied_verifier_sha256':hashlib.sha256(code.read_bytes()).hexdigest(),
             'replay_verification_sha256':hashlib.sha256((out/'verification.json').read_bytes()).hexdigest(),
             'scope':'Numerical replay without Git/network/production imports; no new source outcomes.'}
    (here/'portable_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))


if __name__=='__main__':main()
