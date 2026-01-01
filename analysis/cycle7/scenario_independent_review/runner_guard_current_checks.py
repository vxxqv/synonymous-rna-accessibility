"""Read-only runner guard checks with invented bytes, no landscape inputs."""
import sys
if not __debug__:raise RuntimeError('Optimized reviewer execution refused')
from pathlib import Path
import hashlib
import json
import subprocess
import tempfile


def main():
    here=Path(__file__).resolve().parent
    runner=here.parents[0]/'run_scenarios.py'
    original=runner.read_bytes()
    tests=[]
    with tempfile.TemporaryDirectory(prefix='guard_',dir=here) as directory:
        root=Path(directory)
        snapshot=root/'runner.py';snapshot.write_bytes(original)
        core=root/'scenario_geometry'/'scenario_geometry.py'
        core.parent.mkdir()
        core.write_text("raise RuntimeError('IMPORT_REACHED_UNEXPECTEDLY')\n",encoding='utf-8')
        freeze=root/'freeze.json';result=root/'result.json';out=root/'out.json'
        def command(expected,flag=None):
            cmd=[sys.executable]+([flag] if flag else [])+[str(snapshot),'--freeze',str(freeze),'--result',str(result),'--out',str(out)]
            p=subprocess.run(cmd,text=True,capture_output=True)
            assert p.returncode and expected in p.stderr,(expected,p.stderr)
            assert 'IMPORT_REACHED_UNEXPECTEDLY' not in p.stderr
            tests.append({'test':expected,'optimization':flag,'passed':True})
        # Neither freeze nor result exists for these first guards.
        command('Optimized Python prohibited','-O')
        command('Optimized Python prohibited','-OO')
        out.write_text('retained',encoding='utf-8')
        command('Refusing existing result')
        assert out.read_text(encoding='utf-8')=='retained'
        out.unlink()
        freeze.write_text(json.dumps({'files':{'scenario_geometry/scenario_geometry.py':'0'*64},'primary_result_sha256':'0'*64}),encoding='utf-8')
        command('Changed frozen source')
        freeze.write_text(json.dumps({'files':{'scenario_geometry/scenario_geometry.py':hashlib.sha256(core.read_bytes()).hexdigest()},'primary_result_sha256':'0'*64}),encoding='utf-8')
        result.write_bytes(b'not even JSON; invented guard input')
        command('Changed primary result')
        assert not out.exists()
    assert runner.read_bytes()==original
    receipt={'status':'PASS_FAILURE_GUARDS','runner_sha256':hashlib.sha256(original).hexdigest(),
             'checks':tests,'no_observed_input_read':True,'guard_standin_import_not_reached':True,
             'scope':'Failure guards only; mathematical calculations independently tested in reference_checks.py.',
             'reviewer_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    target=here/'runner_guard_current_receipt.json'
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))


if __name__=='__main__':main()
