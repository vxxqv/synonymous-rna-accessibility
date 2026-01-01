from pathlib import Path
import hashlib,json,math,shutil

root=Path(__file__).resolve().parents[2];base=root/'work/persist_cycle3'
repo=root/'outputs/synonymous-rna-accessibility'
dest=repo/'results/programme/persist_v2/train'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=base/'training_v2/predictions.json';data=json.loads(p.read_text())
assert not data['holdout_targets_read']
assert len(data['development_targets'])==22 and len(data['holdout_ids'])==15 and len(data['predictions'])==6
assert all(len(v)==15 and all(math.isfinite(x) for x in v) for v in data['predictions'].values())
assert data['features_sha256']==sha(base/'features_v2/features.json')
assert data['protocol_sha256']==sha(base/'prediction_v2/protocol.json')
dest.mkdir(parents=True,exist_ok=False)
files={'predictions.json':p,'features.json':base/'features_v2/features.json',
       'training_authorization.json':base/'training_authorization.json',
       'authorize_training.py':base/'authorize_training.py','preserve_training.py':Path(__file__)}
for name,source in files.items():shutil.copyfile(source,dest/name)
receipt={'status':'six_complete_predictions_frozen_before_holdout_access',
 'feature_stage_protocol_commit':'5f321adb2ac37ca519949d8a852ffc0ac4f891b3',
 'sha256':{name:sha(dest/name) for name in files},
 'development_count':22,'holdout_count':15,'holdout_cells_decoded':0,
 'selected_penalties':{k:v['selected_penalty'] for k,v in data['tuning'].items()},
 'scope':'Prediction in reported cellular-coefficient units, one selected reporter context; no physical-rate calibration.'}
(dest/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps(receipt))
