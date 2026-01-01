"""Independent exact arithmetic for the invented source-R collapse fixture."""
import argparse,csv,hashlib,json
from fractions import Fraction as F
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def verify(path,out):
    path=Path(path);out=Path(out)
    if out.exists():raise FileExistsError('Fresh receipt required')
    with path.open(newline='',encoding='utf8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    # Independently collapse the invented two-base allele: median(100,102)=101,
    # median(10,12)=11; the other allele has count30. All use WT101-11-30=60.
    expected={'toy:1:AA:CC':(101,11,60,F(23,121)),'toy:1:A:G':(101,30,60,F(61,121))}
    if len(rows)!=len(expected) or {r['VAR_ID'] for r in rows}!=set(expected):raise AssertionError('Fixture identity loss')
    discrepancies=[]
    for r in rows:
        dp,cao,wt,q=expected[r['VAR_ID']]
        if tuple(int(r[x]) for x in ('DP','CAO','WT'))!=(dp,cao,wt):raise AssertionError('Source collapse disagreement')
        error=abs(F(r['q'])-q)
        if error>F(1,10**14):raise AssertionError('Source floating ratio disagreement')
        discrepancies.append(error)
    receipt={'status':'PASS_INVENTED_SOURCE_NORMALIZATION','real_source_rows_read':0,'checks':10,
             'expected_ratios':{k:str(v[3]) for k,v in expected.items()},
             'maximum_absolute_difference':float(max(discrepancies)),
             'scope':'Source helper numerical collapse versus independent rational invented counts; not historical data equivalence.',
             'input_sha256':sha(path),'verifier_sha256':sha(__file__)}
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8');print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();verify(a.input,a.out)
