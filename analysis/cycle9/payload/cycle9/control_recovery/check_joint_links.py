from pathlib import Path
from itertools import product
from fractions import Fraction as F
import json, hashlib
HERE=Path(__file__).resolve().parent
def main():
 if not __debug__:raise RuntimeError('Optimized Python refused')
 d=json.loads((HERE/'results_v1/result.json').read_text())
 vals={r['id']:F(r['value']['exact']) for s in d['state_results'] for r in s['ratios']}
 cases={r['id']:r for r in d['contrast_quotients']}
 require=lambda ok: None if ok else (_ for _ in ()).throw(ValueError('Joint record disagreement'))
 expected=set(product(*[[f'{s}_R{i}_D{j}' for i,j in product(range(1,4),repeat=2)] for s in ('C','M','E')]))
 got=[]
 for r in d['joint_completions']:
  c,m,e=(r[x+'_ratio_id'] for x in ('C','M','E'));got.append((c,m,e))
  require(F(cases[r['CM_case_id']]['value']['exact'])==vals[c]/vals[m])
  require(F(cases[r['CE_case_id']]['value']['exact'])==vals[c]/vals[e])
  require(r['CM_direction']=='below_one' and r['CE_direction']=='below_one')
 require(len(got)==729 and len(set(got))==729 and set(got)==expected)
 out=HERE/'root_independent_review/joint_links.json';require(not out.exists())
 out.write_text(json.dumps({'status':'PASS','all_joint_links':729,'shared_C_verified':True,'case_quotient_checks':1458,
 'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
 print(out.read_text())
if __name__=='__main__':main()
