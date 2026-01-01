from pathlib import Path
import hashlib,json,math

root=Path(__file__).resolve().parents[2];base=root/'work/persist_cycle3'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
authorization=json.loads((base/'features_authorization.json').read_text())
assert authorization['authorized_stages']==['features']
assert all(sha(root/p)==h for p,h in authorization['input_sha256'].items())
protocol=json.loads((base/'prediction_v2/protocol.json').read_text())
path=base/'features_v2/features.json';features=json.loads(path.read_text())
assert features['protocol_sha256']==sha(base/'prediction_v2/protocol.json')
assert features['ids']==protocol['all_ids'] and len(features['ids'])==37
assert set(features['values'])==set(features['diagnostics'])==set(features['ids'])
assert not list((base/'features_v2').glob('*failure*'))
for key in features['ids']:
    values=features['values'][key];d=features['diagnostics'][key]
    assert len(values)==83 and all(math.isfinite(v) for v in values)
    assert d['rescaling_calls']==1 and not d['clipped'] and d['length']==928
    assert abs(d['value']-d['row_formula_value'])<=1e-12
    assert 0<=d['value']<=1 and values[-1]==d['value']
    assert json.loads((base/'features_v2'/f'{key}_fold.json').read_text())==d
result={**authorization,'authorized_stages':['train'],
 'feature_stage_completed':{'count':37,'feature_sha256':sha(path),'no_failure_or_warning':True},
 'restriction':'Only22 development endpoint cells may be decoded. All15 held-out outcome cells remain sealed until a committed prediction freeze exists.'}
result['input_sha256']={**authorization['input_sha256'],str(path.relative_to(root)):sha(path)}
dest=base/'training_authorization.json';assert not dest.exists()
dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'training_only_authorized','feature_count':37,'feature_sha256':sha(path),
 'MUP_minimum':min(v[-1] for v in features['values'].values()),'MUP_maximum':max(v[-1] for v in features['values'].values()),
 'cellular_target_cells_decoded_before_authorization':0,'authorization_sha256':sha(dest)}))
