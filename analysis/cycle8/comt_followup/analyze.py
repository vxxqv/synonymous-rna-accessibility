"""Exposed COMT diagnostics. Fixed predictions only, no model fitting."""
if not __debug__:raise RuntimeError('Optimized mode refused')
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import argparse,csv,hashlib,json,math
import numpy as np

NAMES=['pooled_median','ROI_median','RNA_PAVA','RNA_linear','RNA_load_linear']
PINS={'normalize/contrasts.json':'5654df32df45840415e483a53127a12744887ac9c6b49a18f0a9387bfa04f083','predict/predictions.json':'a28d2075d078224afd51328056604a145da1f59928486dac8df857c7c2f76bfe','evaluate/evaluation.json':'c3590c1f885323617e402c0b3fb99c8bd879a693aaff47993ef260bc6cb06452'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def exact(v):return Q(v['exact_stored_input_arithmetic'])
def putcsv(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    p=argparse.ArgumentParser();p.add_argument('--actual',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.exists():raise FileExistsError('No overwrite')
    assert all(sha(a.actual/k)==v for k,v in PINS.items())
    contrasts,predictions,evaluation=[json.loads((a.actual/k).read_text()) for k in PINS]
    variants=[];all_scores=[];site_rows=[];influence=[];fits=[];flags=[];opposite={};summary={};maxcond=0.0
    for case in ['code','geo']:
        lookup={r['hgvs_nt']:r for r in contrasts[case]};assert len(lookup)==125
        opposite[case]=[]
        for vid,r in sorted(lookup.items()):
            rr=r['contrasts']['R'];tt=r['contrasts']['T'];disagree=rr['resolved_sign']*tt['resolved_sign']==-1
            margin=lambda z:z['lower'] if z['resolved_sign']==1 else -z['upper'] if z['resolved_sign']==-1 else 0.0
            radius=min(margin(rr),margin(tt))/2 if disagree else 0.0
            row={'map':case,'variant':vid,'site':r['site'],'ROI':r['ROI'],'is_SNV':r['is_SNV'],'opposite_resolved_signs':disagree,'assay_log_widening_radius_exclusive':radius}
            for name in ['R','T','L','L_mid']:
                for field in ['center','lower','upper','resolved_sign']:row[name+'_'+field]=r['contrasts'][name][field]
            variants.append(row)
            if disagree:opposite[case].append(vid)
        summary[case]={}
        for load in ['heavy','mid']:
            ev=evaluation['cases'][case][load];pp=predictions[case][load]
            assert len(pp['predictions'])==125 and not pp['unavailable'] and ev['scores']['sites']==53
            for label,sc in [('all',ev['scores']),*ev['sensitivities'].items()]:
                for model in NAMES:
                    all_scores.append({'map':case,'load':load,'subset':label,'model':model,'variants':sc['variants'],'sites':sc['sites'],'site_MAE_exact_stored_inputs':str(exact(sc['site_MAE'][model])),'site_MAE_nats':float(exact(sc['site_MAE'][model]))})
            errs={site:{m:exact(v[m]) for m in NAMES} for site,v in ev['scores']['per_site_errors'].items()}
            for site,values in sorted(errs.items()):
                members=[r for r in lookup.values() if r['site']==site]
                assert len({r['ROI'] for r in members})==1
                for model in NAMES:site_rows.append({'map':case,'load':load,'site':site,'ROI':members[0]['ROI'],'variants':len(members),'model':model,'MAE_exact_stored_inputs':str(values[model]),'MAE_nats':float(values[model])})
            sums={m:sum(v[m] for v in errs.values()) for m in NAMES};ns=len(errs)
            for m in NAMES:assert sums[m]/ns==exact(ev['scores']['site_MAE'][m])
            diffs=[v['RNA_load_linear']-v['ROI_median'] for v in errs.values()]
            local=[]
            for site,values in sorted(errs.items()):
                scores={m:(sums[m]-values[m])/(ns-1) for m in NAMES}
                nv=125-sum(r['site']==site for r in lookup.values())
                gains={m:1-scores['RNA_load_linear']/scores[m] if scores[m]>0 else None for m in NAMES[:-1]}
                passed=nv>=100 and ns-1>=40 and all(g is not None and g>=Q(1,10) for g in gains.values())
                row={'map':case,'load':load,'omitted_assessment_site':site,'retained_sites':ns-1,'retained_variants':nv,'original_operational_rule_pass':passed}
                for m in NAMES:row[m+'_MAE_nats']=float(scores[m])
                for m,g in gains.items():row['gain_vs_'+m]=None if g is None else float(g)
                influence.append(row);local.append(row)
            flagcounts=Counter()
            for r in pp['predictions']:
                for model,f in r['range_flags'].items():
                    row={'map':case,'load':load,'variant':r['hgvs_nt'],'fold':r['fold'],'model':model,**f}
                    row.setdefault('load_outside',False);flags.append(row)
                    for name,v in f.items():
                        if v:flagcounts[model+'/'+name]+=1
            for fold in pp['folds']:
                assert fold['status']=='predicted';tr=[lookup[v] for v in fold['training_ids']]
                counts=Counter(r['site'] for r in tr);ww=np.array([1/counts[r['site']] for r in tr]);xn=np.array([r['contrasts']['R']['center'] for r in tr]);ln=np.array([r['contrasts']['L' if load=='heavy' else 'L_mid']['center'] for r in tr])
                for model in ['RNA_linear','RNA_load_linear']:
                    f=fold['fit'][model];xx=xn[:,None] if model=='RNA_linear' else np.column_stack([xn,ln]);cc=xx-np.average(xx,axis=0,weights=ww);sv=np.linalg.svd(cc*np.sqrt(ww[:,None]),compute_uv=False)
                    observed=float(sv[0]/sv[-1]) if sv[-1]>0 else None
                    expected=f['condition_number'];difference=abs(observed-expected) if observed is not None and expected is not None else None
                    if difference is not None:
                        assert math.isclose(observed,expected,rel_tol=1e-10,abs_tol=1e-12);maxcond=max(maxcond,difference)
                    fits.append({'map':case,'load':load,'fold':fold['fold'],'model':model,'centered_rank':f['centered_rank'],'design_rank':f['design_rank'],'condition_number_stored':expected,'condition_number_independent_SVD':observed,'condition_difference':difference,'intercept':str(exact(f['intercept'])),'RNA_slope':str(exact(f['rna_slope'])),'load_slope':str(exact(f['load_slope'])),'gram_determinant':str(exact(f['gram_determinant'])),'degeneracy':f['degeneracy']})
            summary[case][load]={'sites_joint_better_than_ROI':sum(d<0 for d in diffs),'sites_joint_worse_than_ROI':sum(d>0 for d in diffs),'sites_equal':sum(d==0 for d in diffs),'single_site_omission_original_rule_passes':sum(r['original_operational_rule_pass'] for r in local),'gain_vs_ROI_omission_min':min(r['gain_vs_ROI_median'] for r in local),'gain_vs_ROI_omission_max':max(r['gain_vs_ROI_median'] for r in local),'range_flags':dict(flagcounts)}
    stable=sorted(set(opposite['code'])&set(opposite['geo']));summary['opposite_resolved_signs']={**opposite,'stable_both_maps':stable};summary['max_absolute_condition_discrepancy']=maxcond
    summary['status']='EXPLORATORY_FIXED_PREDICTION_DIAGNOSTICS_PENDING_SEPARATE_REVIEW'
    summary['limits']='Descriptive conditional proxies; fixed-prediction site omissions are not refits. No confidence, substitution-equivalence, independent replication or kinetic mechanism.'
    assert all(sha(a.actual/k)==v for k,v in PINS.items())
    a.out.mkdir(parents=True)
    for name,rows in [('variant_stability.csv',variants),('all_scores.csv',all_scores),('site_errors.csv',site_rows),('site_omission_influence.csv',influence),('fit_diagnostics.csv',fits),('extrapolation_flags.csv',flags)]:putcsv(a.out/name,rows)
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    (a.out/'manifest.json').write_text(json.dumps({'inputs':PINS,'code_sha256':sha(Path(__file__)),'protocol_sha256':sha(Path(__file__).with_name('protocol.md')),'outputs':{f.name:sha(f) for f in sorted(a.out.iterdir()) if f.is_file()}},indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary))
if __name__=='__main__':main()
