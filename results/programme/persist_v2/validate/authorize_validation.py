"""Authenticate the already committed prediction blob before held-out access."""
from pathlib import Path
import hashlib,json,subprocess

root=Path(__file__).resolve().parents[2];base=root/'work/persist_cycle3'
repo=root/'outputs/synonymous-rna-accessibility'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
commit='a6a2ded62c3fd7581f04cff29d490fbb8796b9d2'
git_path='results/programme/persist_v2/train/predictions.json'
prediction=base/'training_v2/predictions.json'
raw=subprocess.check_output(['git','-C',str(repo),'show',f'{commit}:{git_path}'])
assert hashlib.sha256(raw).hexdigest()==sha(prediction)=='2a22b0977afc031bfd2eceba15f5769daffad51f01a11319e81f787c7090a103'
remote=subprocess.check_output(['git','-C',str(repo),'ls-remote','origin','refs/heads/main']).decode().split()[0]
assert remote==commit
freeze={'prediction_sha256':sha(prediction),'commit':commit,'git_path':git_path,
        'repository':str(repo),'remote_main_verified':remote,
        'holdout_outcomes_decoded_before_this_freeze':0}
file=base/'prediction_freeze.json';assert not file.exists()
file.write_text(json.dumps(freeze,indent=2)+'\n',encoding='utf-8')
auth=json.loads((base/'training_authorization.json').read_text())
assert all(sha(root/p)==h for p,h in auth['input_sha256'].items())
auth.update({'authorized_stages':['validate'],
             'restriction':'Decode exactly15 frozen held-out endpoint cells after authentication of the committed predictions. Retire these outcomes from future untouched validation; no tuning, subset changes or replacement feature.'})
auth['input_sha256']={**auth['input_sha256'],str(file.relative_to(root)):sha(file),str(prediction.relative_to(root)):sha(prediction)}
out=base/'validation_authorization.json';assert not out.exists()
out.write_text(json.dumps(auth,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'frozen_validation_authorized','prediction_commit':commit,'prediction_freeze_sha256':sha(file),
 'authorization_sha256':sha(out),'holdout_cells_decoded':0}))
