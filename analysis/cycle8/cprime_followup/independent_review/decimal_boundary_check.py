"""Exact display-boundary audit, no changes to the prescribed decision rule."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reconstruction',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if not __debug__:raise RuntimeError('Optimized execution prohibited')
    if a.out.exists():raise FileExistsError(a.out)
    pin=sha(a.reconstruction);r=json.loads(a.reconstruction.read_text())
    cases=[]
    for row in r['profile_records']:
        for name in ['A','B']:
            value=F(row[name]);distance=value-2
            if abs(distance)<=F(1,10**12):
                cases.append({'intended':row['intended'],'cell':row['cell'],'quantity':name,
                  'value':str(value),'minus_two_exact':str(distance),
                  'at_least_two':value>=2,'at_most_two':value<=2,
                  'interpretation':'Exact stored-decimal classification only. This closeness cutoff diagnoses display precision and is not a changed biological decision rule.'})
    robust=[x for x in r['decimal_boundary_diagnostics']['transition_and_all_three_margins'] if x['quantity']=='smallest_recorded_assignment_fold_minus_two']
    if len(robust)!=8 or any(F(x['exact_minus_two'])<=F(1,10**12) for x in robust):raise AssertionError('Unexpected all-three threshold dependence')
    result={'status':'PASS','input_sha256':pin,'source_sha256':sha(Path(__file__)),
      'transition_boundary_cases':cases,'point_assignment_cases_near_two':r['decimal_boundary_diagnostics']['cases_with_assignment_folds_near_two'],
      'all_three_robust_cases':8,'smallest_all_three_margin_exact':str(min(F(x['exact_minus_two']) for x in robust)),
      'all_primary_classifications_retained':True,'new_phenotype_evaluation':False}
    if sha(a.reconstruction)!=pin:raise RuntimeError('Input drift')
    a.out.mkdir(parents=True);(a.out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
if __name__=='__main__':main()
