"""Restore immutable audit files and optionally verify COMT in a fresh workspace.

This relocates paths, not scientific bytes. It does not refit the published
classifier or change the frozen prediction procedure. Python 3.12 is used.
"""
import argparse,gzip,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path,PurePosixPath

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def long(p):
    s=str(p.resolve())
    return Path('\\\\?\\'+s) if os.name=='nt' and not s.startswith('\\\\?\\') else Path(s)
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')

def safe_join(root,relative):
    if not isinstance(relative,str) or not relative or '\\' in relative or ':' in relative:
        raise ValueError('Invalid portable relative path')
    path=PurePosixPath(relative)
    if path.is_absolute() or any(x in {'.','..',''} for x in relative.split('/')):
        raise ValueError('Unsafe portable relative path: '+relative)
    joined=root.joinpath(*path.parts)
    if root.resolve() not in joined.resolve().parents:
        raise ValueError('Path escapes intended root')
    for ancestor in [joined,*joined.parents]:
        if ancestor.is_symlink():raise ValueError('Symlinks are not archive inputs')
        if ancestor==root:break
    return joined
def main():
    if not __debug__:raise RuntimeError('Optimized execution is refused')
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--check-comt',action='store_true');a=ap.parse_args()
    repo=long(a.repo);out=long(a.out)
    if out.exists():raise FileExistsError('Fresh output directory required')
    if repo==out or repo in out.parents or out in repo.parents:raise ValueError('Output must be outside repository and its ancestors')
    source=repo/'analysis/cycle8';manifests=sorted(source.glob('VERIFIED_*_CONTENTS.json'))
    if not manifests:raise ValueError('Verified content manifests missing')
    target=out/'work/research_cycle8'
    inventory_path=source/'REPLAY_FILE_INVENTORY.json'
    inventory=json.loads(inventory_path.read_text())
    sealed=[];seen_inventory=set()
    for row in inventory['files']:
        p=safe_join(source,row['path']);dest=safe_join(target,row['path'])
        key=str(dest).casefold() if os.name=='nt' else str(dest)
        if key in seen_inventory:raise ValueError('Duplicate inventory path')
        seen_inventory.add(key)
        if sha(p)!=row['sha256'] or p.stat().st_size!=row['bytes']:raise ValueError('Sealed inventory mismatch: '+row['path'])
        sealed.append((p,dest))
    allentries=[];seen_logical=set()
    for manifest in manifests:
        for e in json.loads(manifest.read_text())['files']:
            if not e['archive_path'].startswith('analysis/cycle8/'):raise ValueError('Unexpected archive root')
            p=safe_join(repo,e['archive_path']);dest=safe_join(target,e['logical_path'])
            key=str(dest).casefold() if os.name=='nt' else str(dest)
            if key in seen_logical:raise ValueError('Duplicate logical destination')
            seen_logical.add(key)
            if e['compression'] not in {'none','gzip'}:raise ValueError('Unknown compression')
            if str(target/p.relative_to(source)).casefold() not in seen_inventory and os.name=='nt':raise ValueError('Logical archive file is not sealed')
            if os.name!='nt' and str(target/p.relative_to(source)) not in seen_inventory:raise ValueError('Logical archive file is not sealed')
            if sha(p)!=e['stored_sha256']:raise ValueError('Archive hash mismatch: '+e['archive_path'])
            allentries.append(e)
    # Only the fully sealed inventory is copied. Bytecode, notebooks and stray
    # local files never enter the reproduction workspace through directory copy.
    target.mkdir(parents=True)
    for original,dest in sealed:
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(original,dest)
    restored=[]
    for e in allentries:
        p=target/e['logical_path'];b=(repo/e['archive_path']).read_bytes()
        if e['compression']=='gzip':b=gzip.decompress(b)
        if hashlib.sha256(b).hexdigest()!=e['sha256'] or len(b)!=e['bytes']:raise ValueError('Logical hash mismatch')
        if p.exists() and p.read_bytes()!=b:raise ValueError('Conflicting logical file')
        p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
        restored.append({'path':e['logical_path'],'sha256':e['sha256']})
    # The original B1 result is stored as a user-readable result in the repo.
    b1=repo/'results/cycle8/cprime/primary_result.json'
    if sha(b1)!='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf':raise ValueError('B1 hash mismatch')
    dest=target/'cprime_results/primary_v1/result.json';dest.parent.mkdir(parents=True);shutil.copyfile(b1,dest)
    # Prior FAS artifacts are needed only by the optional exact certificate check.
    for stored,rel,pin in [('results/cycle6/full_result_v1.json.gz','work/research_cycle6/actual_v1/result.json','3f34f322acc8bc9a6f6b414278d6581c458bce08f0f3bb18a3bcccee7320e466'),('results/cycle7/robustness.json.gz','work/research_cycle7/results/robustness.json','861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa')]:
        b=gzip.decompress((repo/stored).read_bytes())
        if hashlib.sha256(b).hexdigest()!=pin:raise ValueError('Historical FAS hash mismatch')
        dest=out/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b)
    result={'status':'MATERIALIZED_BYTE_VERIFIED','logical_files':len(restored),'sealed_files':len(sealed),'sealed_inventory_sha256':sha(inventory_path),'logical_manifest_hashes':{x.name:sha(x) for x in manifests},'restored':restored,'limits':'Archived historical authorizations retain their original paths. This materializer does not claim new prospective validation or historical RF training reproduction.'}
    if a.check_comt:
        ref=target/'comt_numeric_reference/v3_delta_review/reference_successor/binary_float_successor'
        # Relocate every already authorized role by its original workspace suffix.
        candidates=list(ref.glob('*authorization*.json'))
        candidates=[p for p in candidates if json.loads(p.read_text()).get('schema')=='COMT_independent_binary_float_authorization_v3']
        if len(candidates)!=1:raise ValueError('Ambiguous original COMT authorization')
        original=candidates[0];auth=json.loads(original.read_text());auth['archive_replay_of_authorization_sha256']=sha(original)
        for role,spec in auth['inputs'].items():
            old=spec['path'].replace('\\','/')
            markers=[m for m in ('/work/research_cycle8/','/outputs/synonymous-rna-accessibility/analysis/cycle8/') if m in old]
            if len(markers)!=1:raise ValueError('Unexpected role path: '+role)
            marker=markers[0]
            dest=target/old.split(marker,1)[1]
            if sha(dest)!=spec['sha256']:raise ValueError('Relocated role hash mismatch: '+role)
            spec['path']=str(dest)
        authpath=out/'comt_replay_authorization.json';dump(authpath,auth)
        verifyout=out/'verification/comt';verifyout.parent.mkdir(exist_ok=True)
        log=out/'comt_verification.log'
        with log.open('w',encoding='utf-8') as stream:
            run=subprocess.run([sys.executable,'-B',str(ref/'reference_binary_float.py'),'--authorization',str(authpath),'--out',str(verifyout)],stdout=stream,stderr=subprocess.STDOUT)
        if run.returncode:raise RuntimeError('COMT replay failed; preserved log: '+str(log))
        rr=json.loads((verifyout/'receipt.json').read_text())
        if rr['status']!='PASS_INDEPENDENT_BINARY_FLOAT_CONTRACT_COMPARISON' or rr['primary_operational_expansion_pass'] is not False:raise ValueError('Unexpected scientific status')
        result.update(status='MATERIALIZED_AND_COMT_VERIFIED',comt_checks=rr['checks'],comt_receipt_sha256=sha(verifyout/'receipt.json'))
    dump(out/'materialization_receipt.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in {'restored','limits'}}))
if __name__=='__main__':main()
