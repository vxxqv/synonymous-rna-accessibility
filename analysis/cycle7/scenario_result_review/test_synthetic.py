"""Small direct controls for the separate observed verifier; no data reads."""
import sys
if not __debug__:raise RuntimeError('Optimized reviewer execution refused')
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,subprocess
import verify_observed as v


def main():
    diamond=[{1,2},{0,3},{0,3},{1,2}]
    first=[Q(2),Q(1),Q(0),Q(2)]
    second=[Q(2),Q(0),Q(1),Q(2)]
    assert v.closure(first,diamond)[0][3]==v.closure(second,diamond)[0][3]==1
    assert v.closure([min(a,b) for a,b in zip(first,second)],diamond)[0][3]==0
    assert v.closure([max(a,b) for a,b in zip(first,second)],diamond)[0][3]==1
    series=[{1},{0,2},{1,3},{2}]
    assert v.closure(first,series)[0][3]==v.closure(second,series)[0][3]==0
    assert v.closure([max(a,b) for a,b in zip(first,second)],series)[0][3]==1
    assert v.closure([Q(0),Q(0)],[{1},{0}])[0][1]==0
    assert v.closure([Q(1),Q(1)],[set(),set()])[0][1]==-float('inf')
    a=v.cell('adaptive',Q(0),Q(1),Q(2))
    assert not v.member(a,Q(0)) and v.member(a,Q(1)) and not v.member(a,Q(3,2))
    a=v.cell('fixed',-float('inf'),Q(0),Q(2))
    assert a['nonempty'] and a['width']=='0' and v.member(a,Q(0))
    assert not v.cell('empty',Q(1),Q(1),Q(2))['nonempty']
    flags=[]
    for flag in ['-O','-OO']:
        p=subprocess.run([sys.executable,flag,str(Path(v.__file__).resolve())],text=True,capture_output=True)
        assert p.returncode and 'Optimized execution refused' in p.stderr
        flags.append(flag)
    target=Path(__file__).with_name('synthetic_receipt.json')
    if target.exists():raise FileExistsError(target)
    target.write_text(json.dumps({'status':'PASS_SYNTHETIC_ONLY','controls':['parallel fixed/adaptive gap','series joint/box gap','zero threshold','disconnected graph','closed tie','empty interval','zero-width singleton'],
        'optimization_rejected':flags,'source_sha256':hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest(),'test_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n',encoding='utf-8')
    print('PASS_SYNTHETIC_ONLY')


if __name__=='__main__':main()
