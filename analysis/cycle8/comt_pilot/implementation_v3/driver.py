"""Guarded staged runner. No execution is authorized by merely supplying paths.
Parent must review/freeze an authorization with exact protocol/code/input hashes.
"""
import argparse,csv,hashlib,json,os,subprocess
from pathlib import Path
import selective_reader as reader
import normalize as norm
import prediction

P=Path(__file__).resolve().parent
PHASES=('extract','rf','normalize','predict','evaluate')
LEGACY_PINS=('protocol','amendment','integration_addendum','contract','code_manifest')
COMMON_PINS=LEGACY_PINS+('normalization_amendment','predecessor_code_manifest')
STAGE_FILES={'extract':{'rows_1.tsv','rows_2.tsv','simulation_1.tsv','simulation_2.tsv','truth_1.tsv','truth_2.tsv','receipt.json'},
             'rf':{'row_status_1.tsv','row_status_2.tsv','status.txt'},
             'normalize':{'normalization.json','contrasts.json'},'predict':{'predictions.json'},'evaluate':{'evaluation.json'}}


def dump(path,obj):Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf8')


def extended(path):
    text=str(Path(path).resolve())
    if os.name=='nt' and not text.startswith('\\\\?\\'):
        text='\\\\?\\UNC\\'+text[2:] if text.startswith('\\\\') else '\\\\?\\'+text
    return Path(text)


def runtime_check(pilot,manifest):
    rscript=Path(manifest['Rscript'])
    if reader.sha(rscript)!=manifest['Rscript_sha256']:raise ValueError('R executable changed')
    for name,digest in manifest['package_files'].items():
        if reader.sha(extended(pilot/'runtime/library'/name))!=digest:raise ValueError('Runtime package file changed '+name)
    for name,digest in manifest['helper_source_sha256'].items():
        if reader.sha(pilot/name)!=digest:raise ValueError('Helper changed '+name)
    if reader.sha(pilot/'source_helpers.R')!=manifest['source_helpers_sha256']:raise ValueError('Helper adapter changed')
    return rscript


def seal(out,phase,auth_hash,inputs):
    out=Path(out)
    files={p.name:reader.sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='sealed_manifest.json'}
    payload={'schema':'COMT_A1_stage_seal_v1','phase':phase,'authorization_sha256':auth_hash,'input_hashes':inputs,'files':files}
    dump(out/'sealed_manifest.json',payload)
    return payload


def stage_manifest(path,phase,expected_inputs=None,pin_names=COMMON_PINS):
    m=json.loads(Path(path).read_text())
    if m.get('schema')!='COMT_A1_stage_seal_v1' or m.get('phase')!=phase:raise ValueError('Wrong prior stage')
    if not STAGE_FILES[phase]<=set(m['files']):raise ValueError('Required prior output is not sealed')
    for n,h in m['files'].items():
        if Path(n).name!=n or reader.sha(Path(path).parent/n)!=h:raise ValueError('Sealed prior file changed')
    if expected_inputs is not None:
        for name in pin_names:
            if m['input_hashes'].get(name)!=expected_inputs[name]:raise ValueError('Prior stage provenance mismatch '+name)
    return m


def load_rows(extract_dir):
    rows=[]
    for roi in (1,2):
        with (extract_dir/f'rows_{roi}.tsv').open(encoding='utf8',newline='') as f:rows.extend(csv.DictReader(f,delimiter='\t'))
    return rows


def load_status(rfdir):
    out={}
    for roi in (1,2):
        with (rfdir/f'row_status_{roi}.tsv').open(encoding='utf8',newline='') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                for name in ('frequency_pass','globally_rejected'):
                    if r[name] not in ('TRUE','FALSE'):raise ValueError('Invalid R boolean')
                    r[name]=r[name]=='TRUE'
                if r['record_id'] in out:raise ValueError('Duplicate RF row status')
                out[r['record_id']]=r
    return out


def execute(phase,auth_path,out):
    reader.require_debug()
    out=Path(out).resolve()
    if out.exists():raise FileExistsError('New output directory required before any input decoding')
    auth_path=Path(auth_path);auth=json.loads(auth_path.read_text())
    if auth.get('schema')!='COMT_A1_authorization_v1' or auth.get('phase')!=phase or auth.get('status')!='approved_frozen_normalization_successor':
        raise ValueError('Missing exact frozen phase authorization')
    if auth.get('historical_reproduction_claim') is not False:raise ValueError('Historical-state claim forbidden')
    if auth.get('outcome_access_authorized') is not True:raise ValueError('Outcome access not authorized')
    roles={}
    for role,item in auth['inputs'].items():
        p=Path(item['path']).resolve()
        if reader.sha(p)!=item['sha256']:raise ValueError('Input hash mismatch '+role)
        roles[role]=p
    if phase not in ('normalize','predict','evaluate'):raise ValueError('Successor reuses sealed extraction and RF; no rerun permitted')
    required=set(COMMON_PINS)
    extra={'extract':{'archive','runtime_manifest'},'rf':{'extract_manifest','runtime_manifest'},'normalize':{'extract_manifest','rf_manifest'},'predict':{'normalization_manifest'},'evaluate':{'normalization_manifest','prediction_manifest'}}[phase]
    if not required|extra<=set(roles):raise ValueError('Missing input-role pin')
    code=json.loads(roles['code_manifest'].read_text())
    if code.get('schema')!='COMT_A1_code_freeze_v1':raise ValueError('Code freeze schema')
    for rel,h in code['files'].items():
        p=(P/rel).resolve()
        if reader.sha(p)!=h:raise ValueError('Frozen implementation changed '+rel)
    own={p.name for p in P.glob('*.py')}|{'apply_source_rf.R','../proxy_models.py','../source_helpers.R','../modeling_utils.r','../summarization_utils.r'}
    if not own<=set(code['files']):raise ValueError('Incomplete implementation freeze')
    contract=json.loads(roles['contract'].read_text())
    if contract.get('schema')!='COMT_A1_selective_contract_v1' or contract['expected_variants']!=125 or contract['expected_sites']!=53:raise ValueError('Wrong target contract')
    targets=contract['targets']
    if len(targets)!=125 or len({r['hgvs_nt'] for r in targets})!=125 or len({r['site'] for r in targets})!=53:
        raise ValueError('Fixed target IDs/sites are incomplete')
    if sum(r['is_SNV'] is True for r in targets)!=113 or len({r['site'] for r in targets if r['is_SNV']})!=52:
        raise ValueError('Fixed SNV subset mismatch')
    for name in ('protocol','amendment','integration_addendum','contract','normalization_amendment','predecessor_code_manifest'):
        if code['named_input_hashes'].get(name)!=reader.sha(roles[name]):raise ValueError('Code/contract association changed')
    inputs={k:reader.sha(v) for k,v in roles.items()}
    auth_hash=reader.sha(auth_path)
    if phase=='extract':
        runtime_check(P.parent,json.loads(roles['runtime_manifest'].read_text()))
        reader.extract(roles['archive'],contract,out)
    elif phase=='rf':
        em=stage_manifest(roles['extract_manifest'],'extract',inputs)
        if em['input_hashes'].get('runtime_manifest')!=inputs['runtime_manifest']:raise ValueError('Extraction used another runtime freeze')
        rscript=runtime_check(P.parent,json.loads(roles['runtime_manifest'].read_text()))
        completed=subprocess.run([str(rscript),'--vanilla',str(P/'apply_source_rf.R'),str(P.parent),str(roles['extract_manifest'].parent),str(out)],capture_output=True,text=True,timeout=600)
        out.mkdir(parents=True,exist_ok=True)
        (out/'execution_log.txt').write_text(completed.stdout+'\n'+completed.stderr,encoding='utf8')
        if completed.returncode!=0:raise RuntimeError('RF source reconstruction failed; preserved execution log')
        if (out/'status.txt').read_text().strip()!='PASS_NEW_FIXED_RUNTIME_RF_APPLICATION':raise ValueError('Missing RF completion')
    elif phase=='normalize':
        predecessor=json.loads(roles['predecessor_code_manifest'].read_text())
        oldpins={**predecessor['named_input_hashes'],'code_manifest':inputs['predecessor_code_manifest']}
        for name in LEGACY_PINS[:-1]:
            if oldpins[name]!=inputs[name]:raise ValueError('Successor changed frozen scientific input '+name)
        stage_manifest(roles['extract_manifest'],'extract',oldpins,LEGACY_PINS)
        rm=stage_manifest(roles['rf_manifest'],'rf',oldpins,LEGACY_PINS)
        if rm['input_hashes'].get('extract_manifest')!=inputs['extract_manifest']:raise ValueError('RF used another extraction')
        result=norm.normalized(load_rows(roles['extract_manifest'].parent),load_status(roles['rf_manifest'].parent),contract)
        summaries=norm.contrast_summaries(result,contract)
        out.mkdir(parents=True);dump(out/'normalization.json',result);dump(out/'contrasts.json',summaries)
    elif phase=='predict':
        stage_manifest(roles['normalization_manifest'],'normalize',inputs)
        summaries=json.loads((roles['normalization_manifest'].parent/'contrasts.json').read_text())
        kernel=prediction.load_kernel(P.parent/'proxy_models.py')
        result=prediction.cross_validate(summaries,kernel)
        out.mkdir(parents=True);dump(out/'predictions.json',result)
    elif phase=='evaluate':
        stage_manifest(roles['normalization_manifest'],'normalize',inputs);pm=stage_manifest(roles['prediction_manifest'],'predict',inputs)
        if pm['input_hashes'].get('normalization_manifest')!=reader.sha(roles['normalization_manifest']):raise ValueError('Predictions came from another normalization')
        if auth.get('sealed_prediction_sha256')!=reader.sha(roles['prediction_manifest']):raise ValueError('Prediction seal must be explicitly authorized before scoring')
        summaries=json.loads((roles['normalization_manifest'].parent/'contrasts.json').read_text())
        predictions=json.loads((roles['prediction_manifest'].parent/'predictions.json').read_text())
        kernel=prediction.load_kernel(P.parent/'proxy_models.py')
        result=prediction.evaluate(predictions,summaries,kernel)
        out.mkdir(parents=True);dump(out/'evaluation.json',result)
    return seal(out,phase,auth_hash,inputs)


if __name__=='__main__':
    reader.require_debug()
    p=argparse.ArgumentParser();p.add_argument('--phase',required=True,choices=PHASES);p.add_argument('--authorization',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(execute(a.phase,a.authorization,a.out),indent=2))
