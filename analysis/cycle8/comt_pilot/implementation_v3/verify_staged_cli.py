"""Disposable invented downstream CLI replay; no archive or count input."""
from pathlib import Path
import json,subprocess,sys
import selective_reader as sr
import driver

P=Path(__file__).resolve().parent
OUT=P/'synthetic_stages_v2'


def main():
    if OUT.exists():raise FileExistsError('Fresh synthetic stages required')
    OUT.mkdir();rows=[];targets=[]
    sites=[s for s in range(1,54) for _ in range(2)]+list(range(1,20))
    nonsnv={i for i,s in enumerate(sites) if s==53}
    nonsnv.update(i for i in range(10))
    for i,s in enumerate(sites):
        roi=1 if s<=33 else 2; x=(s%7)-3; z=((i*3)%11)-5;t=2+x+2*z
        cs={n:{'center':float(v),'lower':float(v),'upper':float(v),'complete':True,'resolved_sign':(v>0)-(v<0)} for n,v in [('R',x),('T',t),('L',z),('L_mid',-z)]}
        target={'hgvs_nt':f'INVENTED_{i}','site':str(s),'ROI':roi,'assessment_fold':s%5,'is_SNV':i not in nonsnv}
        targets.append(target);rows.append({**target,'contrasts':cs})
    normdir=OUT/'normalization';normdir.mkdir()
    driver.dump(normdir/'contrasts.json',{'code':rows,'geo':rows})
    driver.dump(normdir/'normalization.json',{'status':'invented fixture, no raw measurement'})
    contract=OUT/'contract.json'
    driver.dump(contract,{'schema':'COMT_A1_selective_contract_v1','expected_variants':125,'expected_sites':53,'targets':targets})
    common={'protocol':P.parent/'protocol_v2.md','amendment':P.parent/'prospective_reconstruction_amendment.md','integration_addendum':P/'integration_addendum.md','contract':contract}
    code={'schema':'COMT_A1_code_freeze_v1','files':{p.name:sr.sha(p) for p in P.glob('*.py')},'named_input_hashes':{k:sr.sha(v) for k,v in common.items()}}
    for rel in ('apply_source_rf.R','../proxy_models.py','../source_helpers.R','../modeling_utils.r','../summarization_utils.r'):code['files'][rel]=sr.sha(P/rel)
    freeze=OUT/'code.json';driver.dump(freeze,code);common['code_manifest']=freeze
    driver.seal(normdir,'normalize','synthetic',{k:sr.sha(v) for k,v in common.items()})
    roles={**common,'normalization_manifest':normdir/'sealed_manifest.json'}
    def auth(phase,paths):
        return {'schema':'COMT_A1_authorization_v1','phase':phase,'status':'approved_frozen_preoutcome','historical_reproduction_claim':False,
                'outcome_access_authorized':True,'synthetic_fixture_only':True,
                'inputs':{k:{'path':str(v.resolve()),'sha256':sr.sha(v)} for k,v in paths.items()}}
    a=OUT/'predict_authorization.json';driver.dump(a,auth('predict',roles))
    command=[sys.executable,str(P/'driver.py'),'--phase','predict','--authorization',str(a),'--out',str(OUT/'prediction')]
    run=subprocess.run(command,capture_output=True,text=True)
    (OUT/'predict_log.txt').write_text(run.stdout+'\n'+run.stderr)
    if run.returncode:raise RuntimeError('Invented prediction CLI failed')
    predseal=OUT/'prediction/sealed_manifest.json';roles['prediction_manifest']=predseal
    e=auth('evaluate',roles);e['sealed_prediction_sha256']=sr.sha(predseal)
    ep=OUT/'evaluate_authorization.json';driver.dump(ep,e)
    command=[sys.executable,str(P/'driver.py'),'--phase','evaluate','--authorization',str(ep),'--out',str(OUT/'evaluation')]
    run=subprocess.run(command,capture_output=True,text=True)
    (OUT/'evaluate_log.txt').write_text(run.stdout+'\n'+run.stderr)
    if run.returncode:raise RuntimeError('Invented evaluation CLI failed')
    result=json.loads((OUT/'evaluation/evaluation.json').read_text())
    if not result['primary_operational_expansion_pass']:raise AssertionError('Known exact load signal should pass')
    for case in ('code','geo'):
        for mode in ('heavy','mid'):
            r=result['cases'][case][mode]
            if r['scores']['variants']!=125 or r['scores']['sites']!=53:raise AssertionError('Coverage loss')
            if r['sensitivities']['SNV']['variants']!=113 or r['sensitivities']['SNV']['sites']!=52:raise AssertionError('SNV denominator')
            if r['scores']['site_MAE']['RNA_load_linear']['float']!=0:raise AssertionError('Signal arithmetic')
    # Explicit seal mismatch is rejected before any scoring result is written.
    e['sealed_prediction_sha256']='0'*64;bad=OUT/'bad_authorization.json';driver.dump(bad,e)
    run=subprocess.run([sys.executable,str(P/'driver.py'),'--phase','evaluate','--authorization',str(bad),'--out',str(OUT/'forbidden')],capture_output=True,text=True)
    if run.returncode==0 or (OUT/'forbidden').exists():raise AssertionError('Bad prediction seal accepted')
    receipt={'status':'PASS_SYNTHETIC_ONLY','real_source_rows_read':0,'targets':125,'sites':53,'SNVs':113,'SNV_sites':52,
      'cases':2,'modes':2,'models':5,'sealed_prediction_then_evaluation':True,'bad_seal_refused_before_output':True,
      'source_sha256':sr.sha(Path(__file__)),'evaluation_sha256':sr.sha(OUT/'evaluation/evaluation.json')}
    driver.dump(OUT/'receipt.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
