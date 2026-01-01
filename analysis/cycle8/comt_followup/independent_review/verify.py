"""Independent fixed-prediction audit. Does not import or execute author analysis."""
if not __debug__:
    raise RuntimeError('Optimized execution refused before input access')
from pathlib import Path
from fractions import Fraction as F
from collections import Counter, defaultdict
import argparse, csv, hashlib, json, math, platform, sys

NAMES=('pooled_median','ROI_median','RNA_PAVA','RNA_linear','RNA_load_linear')
PINS={'normalize/contrasts.json':'5654df32df45840415e483a53127a12744887ac9c6b49a18f0a9387bfa04f083','predict/predictions.json':'a28d2075d078224afd51328056604a145da1f59928486dac8df857c7c2f76bfe','evaluate/evaluation.json':'c3590c1f885323617e402c0b3fb99c8bd879a693aaff47993ef260bc6cb06452'}
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def q(x): return F(x['exact_stored_input_arithmetic'])
def med(v):
    s=sorted(v); n=len(s)
    return s[n//2] if n%2 else (s[n//2-1]+s[n//2])/2
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--actual',type=Path,required=True);ap.add_argument('--followup',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists(): raise FileExistsError(a.out)
    before={str(a.actual/k):sha(a.actual/k) for k in PINS}
    for k,h in PINS.items(): assert before[str(a.actual/k)]==h
    mf=json.loads((a.followup/'manifest.json').read_text())
    for k,h in mf['outputs'].items(): assert sha(a.followup/k)==h
    assert sha(a.followup.parent/'analyze.py')==mf['code_sha256']
    assert sha(a.followup.parent/'protocol.md')==mf['protocol_sha256']
    audit_inputs={**before,**{str(a.followup/k):sha(a.followup/k) for k in mf['outputs']},str(a.followup/'manifest.json'):sha(a.followup/'manifest.json'),str(a.followup.parent/'analyze.py'):mf['code_sha256'],str(a.followup.parent/'protocol.md'):mf['protocol_sha256']}
    c,p,e=[json.loads((a.actual/k).read_text()) for k in PINS]
    csvs={k:list(csv.DictReader((a.followup/(k+'.csv')).open(newline='',encoding='utf-8'))) for k in ('variant_stability','all_scores','site_errors','site_omission_influence','fit_diagnostics','extrapolation_flags')}
    counts=Counter(); maxdiff=defaultdict(float)
    def equal_float(actual,expected,key,tol=2e-12):
        av=float(actual);ev=float(expected);d=abs(av-ev);maxdiff[key]=max(maxdiff[key],d)
        assert math.isclose(av,ev,rel_tol=tol,abs_tol=tol),(key,av,ev)
        counts[key]+=1
    def b(s): assert s in ('True','False');return s=='True'
    vi={(r['map'],r['variant']):r for r in csvs['variant_stability']};assert len(vi)==len(csvs['variant_stability'])==250
    so={(r['map'],r['load'],r['subset'],r['model']):r for r in csvs['all_scores']};assert len(so)==len(csvs['all_scores'])==80
    se={(r['map'],r['load'],r['site'],r['model']):r for r in csvs['site_errors']};assert len(se)==len(csvs['site_errors'])==1060
    oi={(r['map'],r['load'],r['omitted_assessment_site']):r for r in csvs['site_omission_influence']};assert len(oi)==len(csvs['site_omission_influence'])==212
    fi={(r['map'],r['load'],int(r['fold']),r['model']):r for r in csvs['fit_diagnostics']};assert len(fi)==len(csvs['fit_diagnostics'])==40
    fl={(r['map'],r['load'],r['variant'],r['model']):r for r in csvs['extrapolation_flags']};assert len(fl)==len(csvs['extrapolation_flags'])==1500
    summary={}; opposite={}; stable_details=[]; cond_values=[]; omission_counts={}
    for mapping in ('code','geo'):
        rows={r['hgvs_nt']:r for r in c[mapping]};assert len(rows)==125
        assert len({r['site'] for r in rows.values()})==53
        assert sum(r['is_SNV'] for r in rows.values())==113
        summary[mapping]={};opposite[mapping]=[]
        for vid,r in rows.items():
            v=vi[mapping,vid];assert v['site']==r['site'] and int(v['ROI'])==r['ROI'] and b(v['is_SNV'])==r['is_SNV']
            # Independently reduce the already frozen per-preparation log values.
            ar={}
            for name,assay in r['assays'].items():
                mm=[med([F(x) for x in g['values']]) for g in assay['groups']]
                assert len(mm)==4
                ar[name]=(float(med(mm)),float(min(mm)),float(max(mm)))
            pairs={'R':('R','D'),'T':('D','F'),'L':('H','R'),'L_mid':('M','R')}
            signs={};margins={}
            for name,(u,w) in pairs.items():
                zz=r['contrasts'][name];calc=(ar[u][0]-ar[w][0],ar[u][1]-ar[w][2],ar[u][2]-ar[w][1])
                for field,value in zip(('center','lower','upper'),calc):
                    equal_float(zz[field],value,'preparation_contrast_identity')
                    equal_float(v[name+'_'+field],zz[field],'csv_contrast_identity',0)
                low,high=zz['lower'],zz['upper'];signs[name]=1 if low>0 else -1 if high<0 else 0
                assert signs[name]==zz['resolved_sign']==int(v[name+'_resolved_sign'])
                margins[name]=low if signs[name]==1 else -high if signs[name]==-1 else 0
            opp=signs['R']*signs['T']==-1;rad=min(margins['R'],margins['T'])/2 if opp else 0
            assert b(v['opposite_resolved_signs'])==opp
            equal_float(v['assay_log_widening_radius_exclusive'],rad,'radius_identity',0)
            if opp: opposite[mapping].append(vid)
        for load in ('heavy','mid'):
            pp=p[mapping][load]; ee=e['cases'][mapping][load];pred={r['hgvs_nt']:r for r in pp['predictions']}
            assert set(pred)==set(rows) and not pp['unavailable'];assert len(pp['folds'])==5
            errors={vid:{model:abs(q(pred[vid]['predictions'][model])-F(rows[vid]['contrasts']['T']['center'])) for model in NAMES} for vid in rows}
            def agg(ids):
                by=defaultdict(list)
                for vid in ids:by[rows[vid]['site']].append(vid)
                per={s:{m:sum((errors[v][m] for v in vs),F())/len(vs) for m in NAMES} for s,vs in by.items()}
                return per,{m:sum((v[m] for v in per.values()),F())/len(per) for m in NAMES}
            per,scores=agg(rows)
            for subset,sc in [('all',ee['scores']),*ee['sensitivities'].items()]:
                ids=list(rows) if subset=='all' else [v for v,r in rows.items() if r['is_SNV']] if subset=='SNV' else [v for v,r in rows.items() if 'ROI_'+str(r['ROI'])==subset]
                assert ids,(subset,ee['sensitivities'].keys())
                pe,ma=agg(ids);assert len(ids)==sc['variants'] and len(pe)==sc['sites']
                for m in NAMES:
                    z=so[mapping,load,subset,m];assert int(z['variants'])==len(ids) and int(z['sites'])==len(pe)
                    assert ma[m]==q(sc['site_MAE'][m])==F(z['site_MAE_exact_stored_inputs']);equal_float(z['site_MAE_nats'],float(ma[m]),'aggregate_float',0);counts['exact_aggregate']+=1
            for s,vals in per.items():
                vs=[v for v in rows if rows[v]['site']==s]
                for m in NAMES:
                    z=se[mapping,load,s,m];assert vals[m]==F(z['MAE_exact_stored_inputs'])==q(ee['scores']['per_site_errors'][s][m]);assert int(z['variants'])==len(vs)
                    assert int(z['ROI'])==rows[vs[0]]['ROI'];equal_float(z['MAE_nats'],float(vals[m]),'site_float',0);counts['exact_site']+=1
            local=[]
            for s in per:
                remaining=[v for v in rows if rows[v]['site']!=s];_,sc=agg(remaining);g={m:1-sc['RNA_load_linear']/sc[m] if sc[m]>0 else None for m in NAMES[:-1]}
                passed=len(remaining)>=100 and len(per)-1>=40 and all(x is not None and x>=F(1,10) for x in g.values())
                z=oi[mapping,load,s];assert int(z['retained_sites'])==52 and int(z['retained_variants'])==len(remaining) and b(z['original_operational_rule_pass'])==passed
                for m in NAMES:equal_float(z[m+'_MAE_nats'],float(sc[m]),'omission_score',0)
                for m in g:
                    assert g[m] is not None
                    equal_float(z['gain_vs_'+m],float(g[m]),'omission_gain',0)
                local.append((passed,float(g['ROI_median'])));counts['omissions']+=1
            dif=[val['RNA_load_linear']-val['ROI_median'] for val in per.values()]
            flagcounts=Counter()
            for fold in pp['folds']:
                assert fold['status']=='predicted';train=fold['training_ids'];test=[v for v,z in pred.items() if z['fold']==fold['fold']]
                assert not(set(train)&set(test));assert set(train)|set(test)==set(rows)
                assert not({rows[v]['site'] for v in train}&{rows[v]['site'] for v in test})
                loadname='L' if load=='heavy' else 'L_mid';x={v:F(rows[v]['contrasts']['R']['center']) for v in train};y={v:F(rows[v]['contrasts'][loadname]['center']) for v in train}
                for vid in test:
                    assert pred[vid]['fold']==fold['fold'];rv=F(rows[vid]['contrasts']['R']['center']);lv=F(rows[vid]['contrasts'][loadname]['center'])
                    for m in ('RNA_PAVA','RNA_linear','RNA_load_linear'):
                        rf={'rna_outside':rv<min(x.values()) or rv>max(x.values()),'load_outside':m=='RNA_load_linear' and (lv<min(y.values()) or lv>max(y.values()))}
                        z=fl[mapping,load,vid,m];assert int(z['fold'])==fold['fold']
                        for k,value in rf.items():
                            assert b(z[k])==value
                            assert pred[vid]['range_flags'][m].get(k,False)==value
                            if value:flagcounts[m+'/'+k]+=1
                        counts['range_sets']+=1
                ns=Counter(rows[v]['site'] for v in train);w={v:F(1,ns[rows[v]['site']]) for v in train};W=sum(w.values())
                mx=sum(w[v]*x[v] for v in train)/W;my=sum(w[v]*y[v] for v in train)/W
                A=sum(w[v]*(x[v]-mx)**2 for v in train);B=sum(w[v]*(x[v]-mx)*(y[v]-my) for v in train);D=sum(w[v]*(y[v]-my)**2 for v in train)
                for m in ('RNA_linear','RNA_load_linear'):
                    ff=fold['fit'][m];z=fi[mapping,load,fold['fold'],m]
                    det=F() if m=='RNA_linear' else A*D-B*B
                    rank=int(A>0) if m=='RNA_linear' else 2 if det>0 else int(A+D>0)
                    assert rank==ff['centered_rank']==int(z['centered_rank']) and rank+1==ff['design_rank']==int(z['design_rank'])
                    assert det==q(ff['gram_determinant'])==F(z['gram_determinant'])
                    for k in ('intercept','rna_slope','load_slope'): assert q(ff[k])==F(z['RNA_slope' if k=='rna_slope' else k])
                    assert ff['degeneracy']==z['degeneracy']
                    # Independent 2x2 Gram eigenvalue formula, not the author's SVD implementation.
                    if m=='RNA_linear':cond=1.0 if A>0 else None
                    else:
                        af,bf,df=map(float,(A,B,D));largest=(af+df+math.hypot(af-df,2*bf))/2
                        cond=largest/math.sqrt(float(det)) if det>0 else None
                    assert cond is not None and ff['condition_number'] is not None
                    equal_float(ff['condition_number'],cond,'gram_condition',1e-10)
                    equal_float(z['condition_number_stored'],ff['condition_number'],'stored_condition',0)
                    equal_float(z['condition_number_independent_SVD'],cond,'csv_svd_condition',1e-10)
                    equal_float(z['condition_difference'],abs(float(z['condition_number_independent_SVD'])-float(z['condition_number_stored'])),'stored_condition_difference',0)
                    cond_values.append(cond);counts['exact_rank_gram_sets']+=1
            summary[mapping][load]={'sites_joint_better_than_ROI':sum(d<0 for d in dif),'sites_joint_worse_than_ROI':sum(d>0 for d in dif),'sites_equal':sum(d==0 for d in dif),'single_site_omission_original_rule_passes':sum(z[0] for z in local),'gain_vs_ROI_omission_min':min(z[1] for z in local),'gain_vs_ROI_omission_max':max(z[1] for z in local),'range_flags':dict(flagcounts)}
            omission_counts[mapping+'/'+load]=sum(z[0] for z in local)
    for k in opposite:opposite[k].sort()
    stable=sorted(set(opposite['code'])&set(opposite['geo']));summary['opposite_resolved_signs']={**opposite,'stable_both_maps':stable}
    stored=json.loads((a.followup/'summary.json').read_text())
    for mapping in ('code','geo'):assert summary[mapping]==stored[mapping]
    assert summary['opposite_resolved_signs']==stored['opposite_resolved_signs']
    maxrecord=max(float(z['condition_difference']) for z in csvs['fit_diagnostics']);assert maxrecord==stored['max_absolute_condition_discrepancy']
    for vid in stable:stable_details.append({'variant':vid,'site':vi['code',vid]['site'],'radius_code':float(vi['code',vid]['assay_log_widening_radius_exclusive']),'radius_geo':float(vi['geo',vid]['assay_log_widening_radius_exclusive'])})
    # Boundary known answers for unstructured assay-log interval widening.
    for m in (F(1,8),F(1),F(3,2)):
        assert m-2*(m/2)==0 and m-2*(m/2-F(1,1000))>0 and m-2*(m/2+F(1,1000))<0
    for path,h in audit_inputs.items():assert sha(Path(path))==h
    a.out.mkdir(parents=True)
    receipt={'status':'PASS_INDEPENDENT_FIXED_PREDICTION_DIAGNOSTICS','assertion_tested':'All saved diagnostic summaries reconstructed from pinned contrasts and fixed predictions, no refit or normalization rerun','input_hashes_before_and_after':audit_inputs,'code_sha256':sha(Path(__file__)),'argv':sys.argv,'python':platform.python_version(),'counts':dict(counts),'maximum_absolute_numeric_differences':dict(maxdiff),'condition_number_range':[min(cond_values),max(cond_values)],'omission_rule_pass_counts':omission_counts,'stable_opposite_details':stable_details,'stable_opposite_sites':len({r['site'] for r in stable_details}),'csv_rows':{k:len(v) for k,v in csvs.items()},'tolerances':{'exact_fractions_flags_identities':'exact','contrast_float_identity':2e-12,'condition_number_relative_and_absolute':1e-10},'limitations':'Exploratory fixed-prediction diagnostics; no retraining, confidence, causal mechanism, independent preparation validation, assay substitution or calibrated noise.'}
    (a.out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
