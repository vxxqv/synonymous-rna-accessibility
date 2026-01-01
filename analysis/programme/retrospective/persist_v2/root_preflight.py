"""Independent synthetic linear-system check before any PERSIST outcomes."""
from pathlib import Path
import hashlib,importlib.util,json,sys
import numpy as np

if not __debug__:raise RuntimeError('Assertions required')
root=Path(__file__).resolve().parents[2]
source=root/'work/persist_cycle3/prediction_v2'
dest=root/'work/persist_cycle3/root_preflight.json'
assert not dest.exists()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
frozen=json.loads((source/'proposed_freeze.json').read_text())
assert frozen['authorized_stages']==[]
assert all(sha(root/p)==h for p,h in frozen['input_sha256'].items())
spec=importlib.util.spec_from_file_location('persist_candidate',source/'persist_prediction.py')
module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
rng=np.random.default_rng(928037)
maximum=0.;count=0
for n,p in [(5,3),(22,83),(8,2),(8,7)]:
    x=rng.normal(size=(n,p));x[:,0]=3
    if p>3:x[:,-1]=x[:,1]+x[:,2]
    y=rng.normal(size=n)
    test=rng.normal(size=(4,p));test[:,0]=3
    for lam in [.1,1.,10.,100.]:
        model=module.fit_ridge(x,y,lam)
        active=np.ptp(x,axis=0)!=0
        z=(x[:,active]-x[:,active].mean(0))/x[:,active].std(0)
        # Augmented least squares, independent of the implementation's normal equations.
        design=np.vstack([np.column_stack([np.ones(n),z]),
                          np.column_stack([np.zeros(z.shape[1]),np.sqrt(lam)*np.eye(z.shape[1])])])
        target=np.concatenate([y,np.zeros(z.shape[1])])
        coefficients=np.linalg.lstsq(design,target,rcond=None)[0]
        expected=coefficients[0]+((test[:,active]-x[:,active].mean(0))/x[:,active].std(0))@coefficients[1:]
        error=float(np.max(np.abs(expected-module.predict_ridge(model,test))))
        maximum=max(maximum,error);assert error<1e-10
        count+=1
report={'status':'passed','frozen_input_hashes_verified':len(frozen['input_sha256']),
        'independent_augmented_lstsq_cases':count,'maximum_prediction_difference':maximum,
        'real_sequences_folded':0,'real_target_cells_decoded':0,
        'source_sha256':sha(source/'persist_prediction.py'),'root_check_sha256':sha(Path(__file__)),
        'review':'Endpoint, six baselines, fixed selected cohort, group-wise CV, rowwise feature definition, numerical guards and stage gates inspected. Effects remain conditional predictions in deposited units, with no absolute rate calibration or independent-gene claim.'}
dest.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report))
