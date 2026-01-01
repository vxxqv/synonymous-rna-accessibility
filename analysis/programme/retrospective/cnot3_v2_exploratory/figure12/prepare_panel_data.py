"""Extract exact Figure12 data from independently verified frozen JSON only."""
import sys
if not __debug__: raise RuntimeError('Optimization refused before input access')
sys.dont_write_bytecode=True
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np

PINS={'evaluation':'95fb1c15cfe6148c41a543240e1c84db0d54b910cbe974eeec26914202ceb2ac',
      'verification':'9594b34f602488594a855a249a7e34db260de108a7afdf87d1eae3e631840c45',
      'cohort_receipt':'75cef8d4b9954e5f4975e9be077572e6e28217025e56b6d76163b24d12f92318',
      'kinetic_proof':'c9e2a82f62c97949b01f87d541ccdf0ed6ca1e1445f5e9f0cb288471b85cc2b9'}
ORDER=('no_change','baseline_linear','baseline_core','amino_acid_core','codon_core','association','codon_noSE','association_noSE','association_nonlinear')
LABELS=('No change','Baseline linear','Baseline / core','AA / core','Codon / core','Association (primary)','Codon / core, no SE','Association, no SE','Association, nonlinear')
PREDICTIONS='41661565bf1eb49398fa5cc9141811db4b0f81b10ff77f55a075984897718b93'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def ranks(values):
    order=np.argsort(values,kind='stable'); result=np.empty(len(values),dtype=float); i=0
    while i<len(values):
        j=i+1
        while j<len(values) and values[order[j]]==values[order[i]]:j+=1
        result[order[i:j]]=(i+1+j)/2; i=j
    return result

def main():
    p=argparse.ArgumentParser()
    for name in ('evaluation','verification','cohort_receipt','kinetic_proof','renderer','out'):p.add_argument('--'+name.replace('_','-'),required=True)
    args=vars(p.parse_args()); target=Path(args['out'])
    if target.exists():raise FileExistsError(target)
    before={name:sha(args[name]) for name in PINS}
    if before!=PINS:raise ValueError('Frozen numerical/reference source changed')
    sources={name:json.loads(Path(args[name]).read_text(encoding='utf-8')) for name in ('evaluation','verification','cohort_receipt')}
    e,v,c=[sources[name] for name in ('evaluation','verification','cohort_receipt')]
    if v['status']!='passed' or v['failures'] or v['input_sha256']['evaluation']!=PINS['evaluation']:raise ValueError('Independent verification does not cover this evaluation')
    if e['provenance']['predictions_sha256']!=PREDICTIONS or v['input_sha256']['predictions']!=PREDICTIONS:raise ValueError('Prediction identity differs')
    if c['analysis_status']!='post_exposure_exploratory' or c['historical_published_cohort_reproduced'] is not False:raise ValueError('Exploratory status differs')
    n=len(e['rows']); groups=len({row['group_ID'] for row in e['rows']})
    if n!=1112 or groups!=1109 or c['final_eligible_genes']!=5416 or e['primary_success'] is not False:raise ValueError('Frozen counts/conclusion changed')
    models=[{'model':name,'label':label,'MAE':e['metrics'][name]['MAE'],'role':'primary' if name=='association' else 'sensitivity' if i>=6 else 'baseline'} for i,(name,label) in enumerate(zip(ORDER,LABELS))]
    gains=[]
    for model,baseline,label in [('association','codon_core','Primary'),('association_noSE','codon_noSE','No SE'),('association_nonlinear','codon_core','Nonlinear IP')]:
        value=(e['metrics'][baseline]['MAE']-e['metrics'][model]['MAE'])/e['metrics'][baseline]['MAE']
        stored=e['primary_comparisons']['codon_core']['relative_MAE_reduction'] if model=='association' else e['sensitivities'][model]['relative_MAE_reduction']
        if abs(value-stored)>1e-14:raise ValueError('Frozen gain arithmetic differs')
        gains.append({'model':model,'baseline':baseline,'label':label,'relative_gain':stored,'percent_gain':100*stored,
                      'baseline_MAE':e['metrics'][baseline]['MAE'],'model_MAE':e['metrics'][model]['MAE']})
    k=e['kinetic']['entire_eligible']; records=k['rows']
    if len(records)!=5416 or k['n']!=5416 or len({r['GeneID'] for r in records})!=5416:raise ValueError('Full kinetic cohort changed')
    delta=np.array([r['log2_lifetime_change'] for r in records]); contrast=np.array([r['dimensionless_rate_contrast'] for r in records])
    if not np.all(np.isfinite(delta)) or not np.all(np.isfinite(contrast)):raise ValueError('No finite-point dropping is permitted')
    dx,ry=ranks(delta),ranks(contrast); rho=float(np.corrcoef(dx,ry)[0,1])
    frozen_rho=k['spearman']['log2_lifetime_change__dimensionless_rate_contrast']
    if abs(rho-frozen_rho)>1e-12:raise ValueError('Rank display does not reproduce frozen correlation')
    for r in records:
        expected=math.log(2)*k['h_ref']/r['baseline_HL']*(1-2**(-r['log2_lifetime_change']))
        if abs(expected-r['dimensionless_rate_contrast'])>1e-11*max(1,abs(expected)):raise ValueError('Coupling identity differs')
    points=[{'GeneID':r['GeneID'],'log2_lifetime_change':float(x),'dimensionless_rate_contrast':float(y),
             'lifetime_change_rank_percent':float(100*rx/5416),'rate_contrast_rank_percent':float(100*rr/5416)} for r,x,y,rx,rr in zip(records,delta,contrast,dx,ry)]
    examples=[{'setting':'1','b':1,'p':2,'r':0.5,'delta':0},{'setting':'2','b':0,'p':3,'r':0.5,'delta':0.5}]
    for setting in examples:
        b,pth,r,d=(Fraction(str(setting[name])) for name in ('b','p','r','delta'))
        if b+pth!=3 or b+d+r*pth!=2 or min(b,pth,b+d)<0:raise ValueError('Reviewed synthetic counterexample fails')
    result={'schema':'Figure12_frozen_panel_data_v1','source_sha256':before,'renderer_sha256':sha(args['renderer']),
            'preparer_sha256':sha(Path(__file__)),'prediction_sha256':PREDICTIONS,
            'prediction_commit':'b709472c473e1f7ff53a20d8c7fe28ec97b43bf1','historical_Git_reauthenticated_by_renderer':False,
            'analysis_status':'post_exposure_exploratory','primary_success':False,
            'physical_dimensions_mm':[174,164],'minimum_text_pt':8.5,'biological_uncertainty':'Not calibrated; no interval bars',
            'panels':{'a':{'n_holdout_genes':n,'n_holdout_groups':groups,'models':models,'metric':'Equal-gene MAE of natural-log KO half-life; lower is better'},
                      'b':{'gains':gains,'primary_practical_threshold_percent':5,'threshold_is_primary_reference_not_sensitivity_success_claim':True},
                      'c':{'n':5416,'rho':frozen_rho,'h_ref':k['h_ref'],'all_points':points,'rank_ties':'average exact ties; rank divided by5416',
                           'raw_ranges':{'log2_lifetime_change':[float(delta.min()),float(delta.max())],'dimensionless_rate_contrast':[float(contrast.min()),float(contrast.max())]},
                           'coupling':'R = ln(2) * h_ref / t0 * (1 - 2**(-Delta)); Delta=log2(t1/t0)','points_excluded':0},
                      'd':{'synthetic_only':True,'k0':3,'k1':2,'r':0.5,'units':'arbitrary consistent rate units','examples':examples,'claim':'Distinct p values with identical total rates; no fitted pathway attribution'}}}
    if {name:sha(args[name]) for name in PINS}!=before:raise ValueError('Source inputs changed')
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('x',encoding='utf-8',newline='\n') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({'panel_data_sha256':sha(target),'all_points':len(points),'models':len(models),'rho':frozen_rho}))

if __name__=='__main__':main()
