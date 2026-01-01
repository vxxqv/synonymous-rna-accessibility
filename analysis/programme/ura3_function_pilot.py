"""Frozen descriptive URA3 square benchmark. No measurement-error inference."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from measured_data_audit import xlsx_rows, dump, sha, write_csv

GRID=[-1000,-300,-100,-30,-10,-3,-1,-.3,0,.3,1,3,10,30,100,300,1000]

def boxcox_prediction(a,b,c):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    if c==0:return a+b
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        v=np.expm1(c*a)+np.expm1(c*b)+1
        return np.where((v>0)&np.isfinite(v),np.log(v)/c,np.nan)

def obstruction(y):
    wt,a,b,ab=map(float,y)
    return max(0.,min(wt,ab)-max(a,b),min(a,b)-max(wt,ab))

def isotonic_fit(x,y):
    ux,inv=np.unique(x,return_inverse=True);total=np.bincount(inv,weights=y);weight=np.bincount(inv)
    blocks=[]
    for i,(s,n) in enumerate(zip(total,weight)):
        blocks.append([i,i,float(s),int(n)])
        while len(blocks)>1 and blocks[-2][2]/blocks[-2][3]>blocks[-1][2]/blocks[-1][3]:
            b=blocks.pop();a=blocks.pop();blocks.append([a[0],b[1],a[2]+b[2],a[3]+b[3]])
    fit=np.zeros(len(ux))
    for lo,hi,s,n in blocks:fit[lo:hi+1]=s/n
    return ux,fit,len(blocks)

def features(row,a,b):
    p,q=int(row['position_A']),int(row['position_B']);L=int(row['CDS_nt'])
    muts=[row['single_A'],row['single_B']]
    gc=sum((m[-1] in 'GC')-(m[0] in 'GC') for m in muts)
    transition=sum({m[0],m[-1]} in [{'A','G'},{'C','T'}] for m in muts)
    return [a+b,(p+q)/(2*L),abs(p-q)/L,gc,transition,float(row['same_codon']=='True'),float(row['stop_codon_changed']=='True')]

def ridge_fit(x,y):
    mean=x.mean(axis=0);scale=x.std(axis=0);scale=np.where(scale>1e-14,scale,1.)
    z=np.column_stack([np.ones(len(x)),(x-mean)/scale]);pen=np.eye(z.shape[1]);pen[0,0]=0
    coef=np.linalg.solve(z.T@z+pen,z.T@y)
    return mean,scale,coef

def controls():
    additive=[0.,.02,-.01,.01]
    assert abs(additive[3]-additive[1]-additive[2])<1e-14 and obstruction(additive)==0
    latent_a=np.array([.01,.02,.04]);latent_b=np.array([.02,.03,.01]);c=3
    a=np.log1p(c*latent_a)/c;b=np.log1p(c*latent_b)/c;ab=np.log1p(c*(latent_a+latent_b))/c
    assert np.max(np.abs(boxcox_prediction(a,b,c)-ab))<1e-14
    assert np.all(np.abs(ab-a-b)>1e-5)
    assert all(obstruction([0,x,y,z])==0 for x,y,z in zip(a,b,ab))
    positive=[0.,math.log(.8),math.log(.8),0.]
    assert abs(obstruction(positive)+math.log(.8))<1e-14
    # A nonmonotone global function 1-.2*z*(2-z), z=x1+x2, gives this same square.
    assert [1-.2*z*(2-z) for z in [0,1,1,2]]==[1.,.8,.8,1.]
    for v in np.linspace(-.1,.1,11):
        for w in np.linspace(-.1,.1,11):
            assert obstruction([0,v,w,v+w])<1e-14
    ux,fit,n=isotonic_fit(np.array([0.,0.,1.,2.]),np.array([2.,0.,0.,3.]))
    assert np.allclose(fit,[2/3,2/3,3]) and n==2
    return {'additive_negative_passed':True,'nonlinear_monotone_additive_false_interaction_passed':True,
            'known_square_obstruction_passed':True,'nonmonotone_global_counterexample_passed':True,
            'isotonic_known_solution_passed':True,'measurement_uncertainty_recovered':False}

def metrics(actual,pred,rows):
    p=np.asarray(pred);good=np.isfinite(p);n=len(actual)
    result={'n':n,'valid_predictions':int(good.sum()),'coverage':float(good.mean())}
    if not good.all():return {**result,'log_MAE':None,'log_RMSE':None,'failure':'Prediction-domain failure; incomplete comparisons are not called superior.'}
    err=p-actual;ae=np.abs(err);by_node=defaultdict(list);by_band=defaultdict(list);by_pair=defaultdict(list)
    for i,r in enumerate(rows):
        for m in [r['single_A'],r['single_B']]:by_node[m].append(float(ae[i]))
        by_band[(min(int(r['position_A']),int(r['position_B']))-1)//50].append(float(ae[i]))
        by_pair[tuple(sorted([int(r['position_A']),int(r['position_B'])]))].append(float(ae[i]))
    return {**result,'log_MAE':float(ae.mean()),'log_RMSE':float(np.sqrt(np.mean(err**2))),
            'single_mutant_macro_MAE':float(np.mean([np.mean(v) for v in by_node.values()])),
            'position_pair_macro_MAE':float(np.mean([np.mean(v) for v in by_pair.values()])),
            'band_macro_MAE':float(np.mean([np.mean(v) for v in by_band.values()])),
            'band_MAE':{str(k):float(np.mean(v)) for k,v in sorted(by_band.items())}}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--work',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);dump(a.out/'synthetic_controls.json',controls())
    repo=Path(__file__).resolve().parents[2]
    source=a.work/'programme_sources/gfp_supplementary_data.xlsx'
    assert sha(source)=='7b612a56351fd4815397df67325fef5e94e7adf3910303a03114dc9aacb4a1b6'
    mapping=list(csv.DictReader((a.audit/'synonymous_group_mapping.csv').open(encoding='utf-8',newline='')))
    mapping=[r for r in mapping if r['column']=='J']
    dev=[r for r in mapping if r['planned_partition_J_only']=='development'];val=[r for r in mapping if r['planned_partition_J_only']=='validation']
    assert len(dev)==340 and len(val)==203 and len(mapping)==684
    assert not ({r[k] for r in dev for k in ['single_A','single_B']}&{r[k] for r in val for k in ['single_A','single_B']})
    cells={rn:d for rn,d in xlsx_rows(source,7)};access=[];exclusions=[]
    def read(rn,phase):
        rn=int(rn);cell=cells[rn].get('J');access.append({'cell':'J'+str(rn),'phase':phase})
        if cell is None or cell['kind']!='n':raise ValueError('Missing numeric fitness '+str(rn))
        x=float(cell['text'])
        if not math.isfinite(x) or x<=0:raise ValueError('Nonpositive/nonfinite reported fitness '+str(rn))
        return math.log(x)
    def inputs(rows,phase):return np.array([[read(r['row_A'],phase),read(r['row_B'],phase)] for r in rows])
    ds=inputs(dev,'development_single');dy=np.array([read(r['row_AB'],'development_double') for r in dev])
    # Only validation single measurements are accessed before the prediction seal.
    vs=inputs(val,'validation_single_input')
    grid=[]
    for c in GRID:
        pred=boxcox_prediction(ds[:,0],ds[:,1],c)
        loss=float(np.mean((pred-dy)**2)) if np.isfinite(pred).all() else None
        grid.append({'c':c,'development_log_MSE':loss})
    best=min([r for r in grid if r['development_log_MSE'] is not None],key=lambda r:(r['development_log_MSE'],abs(r['c']),r['c']))['c']
    ux,iso,plateaus=isotonic_fit(ds.sum(axis=1),dy)
    dx=np.array([features(r,*s) for r,s in zip(dev,ds)]);vx=np.array([features(r,*s) for r,s in zip(val,vs)])
    mean,scale,coef=ridge_fit(dx,dy)
    def predict(s,x):
        with np.errstate(invalid='ignore',divide='ignore'):
            linear_w=np.log(np.exp(s[:,0])+np.exp(s[:,1])-1)
        return {'additive_log':s.sum(axis=1),'additive_fitness':linear_w,
                'global_BoxCox':boxcox_prediction(s[:,0],s[:,1],best),
                'isotonic_one_index':np.interp(s.sum(axis=1),ux,iso),
                'sequence_position_ridge':np.column_stack([np.ones(len(s)),(x-mean)/scale])@coef}
    predictions=predict(vs,vx);devpred=predict(ds,dx)
    safe=lambda x:float(x) if math.isfinite(float(x)) else None
    seal={'scope':'Predictions written before validation double values are accessed','source_sha256':sha(source),
          'protocol_sha256':sha(repo/'analysis/programme/measured_function_protocol_v1.json'),
          'eligibility_sha256':sha(repo/'analysis/programme/measured_function_eligibility_v3.json'),
          'model_c':best,'isotonic_plateaus':plateaus,'development_grid':grid,
          'validation':[{'group':r['group'],'row_AB':int(r['row_AB']),'single_log_fitness':list(map(float,s)),
                         'predictions':{k:safe(v[i]) for k,v in predictions.items()}} for i,(r,s) in enumerate(zip(val,vs))]}
    seal_path=a.out/'sealed_validation_predictions.json';dump(seal_path,seal);seal_hash=sha(seal_path)
    # No model, threshold or exclusion choice follows from the validation targets.
    vy=np.array([read(r['row_AB'],'validation_double_after_prediction_seal') for r in val])
    output=[];summaries={}
    for partition,rows,s,y,pred in [('development',dev,ds,dy,devpred),('validation',val,vs,vy,predictions)]:
        margins=[];contrasts=[];fit_contrasts=[]
        for i,(r,ss,yy) in enumerate(zip(rows,s,y)):
            B=obstruction([0,*ss,yy]);margins.append(B);contrasts.append(yy-ss.sum())
            fc=math.exp(yy)-math.exp(ss[0])-math.exp(ss[1])+1;fit_contrasts.append(fc)
            output.append({'partition':partition,'group':r['group'],'single_A':r['single_A'],'single_B':r['single_B'],
                 'position_A':r['position_A'],'position_B':r['position_B'],'w_WT':1.,'w_A':math.exp(ss[0]),'w_B':math.exp(ss[1]),'w_AB':math.exp(yy),
                 'log_interaction':float(yy-ss.sum()),'fitness_interaction':fc,'square_obstruction_log':B,
                 'conservative_uniform_four_cell_log_error_radius':B/2,
                 **{k+'_prediction_log':safe(v[i]) for k,v in pred.items()}})
        margins=np.asarray(margins);contrasts=np.asarray(contrasts);fc=np.asarray(fit_contrasts)
        summaries[partition]={'groups':len(rows),'unique_singles':len({r[k] for r in rows for k in ['single_A','single_B']}),
             'unique_position_pairs':len({tuple(sorted([int(r['position_A']),int(r['position_B'])])) for r in rows}),
             'log_interaction_mean':float(contrasts.mean()),'log_interaction_median':float(np.median(contrasts)),
             'fitness_interaction_mean':float(fc.mean()),'fitness_interaction_median':float(np.median(fc)),
             'reported_w_below_source_cutoff_0p5':int(sum(row[k]<.5 for row in output if row['partition']==partition for k in ['w_A','w_B','w_AB'])),
             'strict_obstruction_counts':{str(tol):int((margins>tol).sum()) for tol in [0,1e-12,1e-8]},
             'box_robust_obstruction_counts':{str(eps):int((margins>2*eps).sum()) for eps in [0,.0001,.001,.005,.01,.025,.05,.1]},
             'maximum_obstruction':float(margins.max()),'prediction_errors':{k:metrics(y,v,rows) for k,v in pred.items()}}
        strata={}
        for label,flag in [('within_codon',True),('between_codons',False)]:
            ix=[i for i,r in enumerate(rows) if (r['same_codon']=='True')==flag]
            strata[label]={'groups':len(ix),'log_interaction_mean':float(contrasts[ix].mean()) if ix else None,
                          'obstructions_above_1e_8':int((margins[ix]>1e-8).sum())}
        summaries[partition]['codon_strata_descriptive']=strata
    write_csv(a.out/'measured_square_results.csv',output);dump(a.out/'cell_access_log.json',access)
    dump(a.out/'results.json',{'source_sha256':sha(source),'prediction_seal_sha256':seal_hash,'selected_BoxCox_c':best,
         'isotonic_plateaus':plateaus,'ridge_lambda':1,'summaries':summaries,
         'numerically_untouched':'All other assay columns; 141 quarantined double outcomes; mitochondrial functional measurements',
         'uncertainty':'Descriptive and deterministic sensitivity only. No calibrated measurement error or significance; shared WT, reporter, selected library and within-set mutants remain dependent.',
         'mechanism':'Neither numerical-scale epistasis nor a four-state obstruction identifies RNA folding; nonmonotone global responses and experimental artifacts remain alternatives.'})
    print(json.dumps({'selected_c':best,'prediction_seal_sha256':seal_hash,'summaries':summaries},allow_nan=False))
if __name__=='__main__':main()
