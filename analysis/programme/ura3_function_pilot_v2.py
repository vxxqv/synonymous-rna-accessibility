"""Honor frozen finite-positive eligibility; retain original failed implementation."""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from measured_data_audit import xlsx_rows, dump, sha, write_csv
from ura3_function_pilot import GRID, boxcox_prediction, obstruction, isotonic_fit, features, ridge_fit, controls, metrics

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--work',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);dump(a.out/'synthetic_controls.json',controls())
    repo=Path(__file__).resolve().parents[2];source=a.work/'programme_sources/gfp_supplementary_data.xlsx'
    assert sha(source)=='7b612a56351fd4815397df67325fef5e94e7adf3910303a03114dc9aacb4a1b6'
    allrows=list(csv.DictReader((a.audit/'synonymous_group_mapping.csv').open(encoding='utf-8',newline='')))
    rows=[r for r in allrows if r['column']=='J'];dc=[r for r in rows if r['planned_partition_J_only']=='development'];vc=[r for r in rows if r['planned_partition_J_only']=='validation']
    assert (len(rows),len(dc),len(vc))==(684,340,203)
    assert not ({r[k] for r in dc for k in ['single_A','single_B']}&{r[k] for r in vc for k in ['single_A','single_B']})
    cells=dict(xlsx_rows(source,7));access=[];exclusions=[]
    def read(rn,phase):
        cell=cells[int(rn)].get('J');access.append({'cell':'J'+str(rn),'phase':phase})
        if cell is None or cell['kind']!='n':return None
        try:x=float(cell['text'])
        except ValueError:return None
        if not math.isfinite(x) or x<=0:return None
        if x<.5:raise ValueError('Positive fitness below documented cutoff requires source reconciliation: J'+str(rn))
        return math.log(x)
    def prepare(candidate,part,read_double):
        retained=[];ss=[];ys=[]
        for r in candidate:
            vals=[read(r['row_A'],part+'_single_A'),read(r['row_B'],part+'_single_B')]
            if read_double:vals.append(read(r['row_AB'],part+'_double'))
            if any(v is None for v in vals):
                exclusions.append({'partition':part,'group':r['group'],'reason':'Nonpositive/nonfinite/unusable reported measurement; no biological meaning imputed',
                                   'unusable_cells':','.join('J'+r[k] for k,v in zip(['row_A','row_B','row_AB'],vals) if v is None)})
            else:
                retained.append(r);ss.append(vals[:2])
                if read_double:ys.append(vals[2])
        return retained,np.asarray(ss),np.asarray(ys)
    dev,ds,dy=prepare(dc,'development',True);vinput,vs,_=prepare(vc,'validation_input',False)
    if len(dev)<100 or len(vinput)<50:
        dump(a.out/'eligibility_failure.json',{'development':len(dev),'validation_inputs':len(vinput),'exclusions':exclusions});raise ValueError('Frozen minimum coverage failed before prediction')
    grid=[]
    for c in GRID:
        pred=boxcox_prediction(ds[:,0],ds[:,1],c);loss=float(np.mean((pred-dy)**2)) if np.isfinite(pred).all() else None
        grid.append({'c':c,'development_log_MSE':loss})
    best=min([r for r in grid if r['development_log_MSE'] is not None],key=lambda r:(r['development_log_MSE'],abs(r['c']),r['c']))['c']
    ux,iso,plateaus=isotonic_fit(ds.sum(axis=1),dy)
    dx=np.array([features(r,*s) for r,s in zip(dev,ds)]);vx=np.array([features(r,*s) for r,s in zip(vinput,vs)]);mean,scale,coef=ridge_fit(dx,dy)
    def predict(s,x):
        with np.errstate(invalid='ignore',divide='ignore'):linear=np.log(np.exp(s[:,0])+np.exp(s[:,1])-1)
        return {'additive_log':s.sum(axis=1),'additive_fitness':linear,'global_BoxCox':boxcox_prediction(s[:,0],s[:,1],best),
                'isotonic_one_index':np.interp(s.sum(axis=1),ux,iso),
                'sequence_position_ridge':np.column_stack([np.ones(len(s)),(x-mean)/scale])@coef}
    dp=predict(ds,dx);vp=predict(vs,vx);safe=lambda x:float(x) if math.isfinite(float(x)) else None
    seal={'scope':'Written before validation double outcome extraction; positive-finite eligibility fixed earlier',
          'source_sha256':sha(source),'endpoint_decision_sha256':sha(repo/'analysis/programme/endpoint_decision_v5.json'),
          'BoxCox_c':best,'grid':grid,'isotonic_plateaus':plateaus,'ridge_lambda':1,
          'validation':[{'group':r['group'],'row_AB':int(r['row_AB']),'single_log_fitness':list(map(float,s)),
                         'predictions':{k:safe(v[i]) for k,v in vp.items()}} for i,(r,s) in enumerate(zip(vinput,vs))]}
    sealpath=a.out/'sealed_validation_predictions.json';dump(sealpath,seal);sealhash=sha(sealpath)
    kept=[];vy=[]
    for i,r in enumerate(vinput):
        y=read(r['row_AB'],'validation_double_after_prediction_seal')
        if y is None:exclusions.append({'partition':'validation_outcome','group':r['group'],'reason':'Unusable reported double outcome; no zero or missing-path inference','unusable_cells':'J'+r['row_AB']})
        else:kept.append(i);vy.append(y)
    val=[vinput[i] for i in kept];vs=vs[kept];vy=np.asarray(vy);vp={k:v[kept] for k,v in vp.items()}
    if len(val)<50:
        dump(a.out/'eligibility_failure.json',{'validation_complete':len(val),'prediction_seal_sha256':sealhash,'exclusions':exclusions});raise ValueError('Frozen validation coverage failed after endpoint eligibility')
    output=[];summaries={}
    for part,rr,s,y,pred in [('development',dev,ds,dy,dp),('validation',val,vs,vy,vp)]:
        margins=[];ee=[]
        for i,(r,ss,yy) in enumerate(zip(rr,s,y)):
            B=obstruction([0,*ss,yy]);margins.append(B);ee.append(yy-ss.sum())
            output.append({'partition':part,'group':r['group'],'single_A':r['single_A'],'single_B':r['single_B'],
                           'position_A':r['position_A'],'position_B':r['position_B'],'same_codon':r['same_codon'],
                           'w_WT':1.,'w_A':math.exp(ss[0]),'w_B':math.exp(ss[1]),'w_AB':math.exp(yy),
                           'log_interaction':float(yy-ss.sum()),'fitness_interaction':math.exp(yy)-math.exp(ss[0])-math.exp(ss[1])+1,
                           'square_obstruction_log':B,'conservative_uniform_four_cell_log_error_radius':B/2,
                           **{k+'_prediction_log':safe(v[i]) for k,v in pred.items()}})
        margins=np.asarray(margins);ee=np.asarray(ee)
        strata={}
        for label,flag in [('within_codon',True),('between_codons',False)]:
            ix=[i for i,r in enumerate(rr) if (r['same_codon']=='True')==flag]
            strata[label]={'n':len(ix),'mean_log_interaction':float(ee[ix].mean()) if ix else None,'obstructions_above_1e_8':int((margins[ix]>1e-8).sum())}
        summaries[part]={'candidate_groups':len(dc if part=='development' else vc),'usable_groups':len(rr),
                        'unique_singles':len({r[k] for r in rr for k in ['single_A','single_B']}),
                        'unique_position_pairs':len({tuple(sorted([int(r['position_A']),int(r['position_B'])])) for r in rr}),
                        'log_interaction_mean':float(ee.mean()),'log_interaction_median':float(np.median(ee)),
                        'strict_obstructions':{str(t):int((margins>t).sum()) for t in [0,1e-12,1e-8]},
                        'box_robust_obstructions':{str(e):int((margins>2*e).sum()) for e in [0,.0001,.001,.005,.01,.025,.05,.1]},
                        'maximum_obstruction':float(margins.max()),'codon_strata':strata,
                        'prediction_errors':{k:metrics(y,v,rr) for k,v in pred.items()}}
    write_csv(a.out/'measured_square_results.csv',output)
    if exclusions:write_csv(a.out/'excluded_groups.csv',exclusions)
    dump(a.out/'cell_access_log.json',access)
    dump(a.out/'results.json',{'prediction_seal_sha256':sealhash,'source_sha256':sha(source),'BoxCox_c':best,'isotonic_plateaus':plateaus,
                             'development_grid':grid,'summaries':summaries,'exclusion_records':len(exclusions),
                             'numerically_untouched':'141 quarantined doubles; all other assay columns except the published C315T scale check; all mitochondrial functional measurements',
                             'scope':'Conditional descriptive point-estimate benchmark. No measurement confidence intervals, independent biological validation, RNA-mechanism claim or full-graph connectivity claim.'})
    print(json.dumps({'BoxCox_c':best,'summaries':summaries,'exclusions':len(exclusions)},allow_nan=False))
if __name__=='__main__':main()
