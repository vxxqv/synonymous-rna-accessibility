"""Frozen five-fold proxy prediction/scoring. No raw-source access."""
from fractions import Fraction as F
import hashlib,importlib.util

KERNEL_SHA='8811af9750e813c798122c7e936b6279a0f3b0e53fea1924cb2d633162a34093'
NAMES=('pooled_median','ROI_median','RNA_PAVA','RNA_linear','RNA_load_linear')


def load_kernel(path):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=KERNEL_SHA:raise ValueError('Held model kernel changed')
    s=importlib.util.spec_from_file_location('comt_held_kernel',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def serial(value):
    if isinstance(value,F):return {'exact_stored_input_arithmetic':str(value),'float':float(value)}
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    return value


def fit_predict(train,assessment,load_name,kernel):
    """Assessment input contains no target field. Fitting accepts training only."""
    if not train:raise ValueError('No complete training examples')
    if any(set(r)!={'hgvs_nt','site','ROI','R','L'} for r in assessment):raise ValueError('Assessment schema must exclude target')
    w=kernel.site_weights([r['site'] for r in train]);x=[r['R'] for r in train];z=[r['L'] for r in train];y=[r['T'] for r in train]
    fit={'pooled_median':kernel.weighted_median(y,w),'ROI_median':kernel.fit_roi_medians(y,w,[r['ROI'] for r in train]),
         'RNA_PAVA':kernel.fit_pava(x,y,w),'RNA_linear':kernel.fit_linear(x,y,w),
         'RNA_load_linear':kernel.fit_linear(x,y,w,z)}
    xx=[r['R'] for r in assessment];zz=[r['L'] for r in assessment]
    pred={'pooled_median':[fit['pooled_median']]*len(xx),'ROI_median':kernel.predict_roi_medians(fit['ROI_median'],[r['ROI'] for r in assessment]),
          'RNA_PAVA':kernel.predict_pava(fit['RNA_PAVA'],xx),'RNA_linear':kernel.predict_linear(fit['RNA_linear'],xx),
          'RNA_load_linear':kernel.predict_linear(fit['RNA_load_linear'],xx,zz)}
    flags={'RNA_PAVA':kernel.range_flags(fit['RNA_PAVA'],xx),'RNA_linear':kernel.range_flags(fit['RNA_linear'],xx),
           'RNA_load_linear':kernel.range_flags(fit['RNA_load_linear'],xx,zz)}
    rows=[]
    for i,r in enumerate(assessment):
        rows.append({**r,'predictions':{n:pred[n][i] for n in NAMES},'range_flags':{n:f[i] for n,f in flags.items()}})
    return {'fit':fit,'predictions':rows,'load':load_name,'training_sites':sorted({r['site'] for r in train}),
            'training_ids':[r['hgvs_nt'] for r in train],'weights':w}


def cross_validate(summaries,kernel):
    result={}
    for case,rows in summaries.items():
        result[case]={}
        for mode,load in [('heavy','L'),('mid','L_mid')]:
            eligible=[r for r in rows if all(r['contrasts'][k]['complete'] for k in ('R','T',load))]
            unavailable=[{'hgvs_nt':r['hgvs_nt'],'site':r['site'],'ROI':r['ROI'],
                          'reason':'incomplete_required_assay_preparations',
                          'missing_contrasts':[k for k in ('R','T',load) if not r['contrasts'][k]['complete']]} for r in rows if r not in eligible]
            folds=[];preds=[]
            for fold in range(5):
                tr=[r for r in eligible if r['assessment_fold']!=fold];te=[r for r in eligible if r['assessment_fold']==fold]
                if {r['site'] for r in tr}&{r['site'] for r in te}:raise ValueError('Codon-site leakage')
                train=[{'hgvs_nt':r['hgvs_nt'],'site':r['site'],'ROI':'ROI_'+str(r['ROI']),
                        'R':r['contrasts']['R']['center'],'L':r['contrasts'][load]['center'],'T':r['contrasts']['T']['center']} for r in tr]
                assessment=[{'hgvs_nt':r['hgvs_nt'],'site':r['site'],'ROI':'ROI_'+str(r['ROI']),
                             'R':r['contrasts']['R']['center'],'L':r['contrasts'][load]['center']} for r in te]
                if not assessment:
                    folds.append({'fold':fold,'status':'no_complete_assessment_rows','training_count':len(train),'assessment_count':0});continue
                try:fitted=fit_predict(train,assessment,load,kernel)
                except ValueError as e:
                    folds.append({'fold':fold,'status':'failed_training_design','reason':str(e),'training_count':len(train),'assessment_count':len(assessment)})
                    unavailable.extend({'hgvs_nt':r['hgvs_nt'],'site':r['site'],'ROI':r['ROI'],'reason':'failed_training_design:'+str(e)} for r in te);continue
                for p in fitted.pop('predictions'):p['fold']=fold;preds.append(p)
                folds.append({'fold':fold,'status':'predicted','training_count':len(train),'assessment_count':len(assessment),**fitted})
            if len({p['hgvs_nt'] for p in preds})!=len(preds):raise ValueError('Duplicate assessment prediction')
            if len(preds)+len(unavailable)!=len(rows):raise ValueError('Lost fixed target')
            result[case][mode]={'folds':folds,'predictions':preds,'unavailable':unavailable,'fixed_variant_denominator':len(rows),'load':load}
    return serial(result)


def exact_prediction(p,name):return F(p['predictions'][name]['exact_stored_input_arithmetic'])


def evaluate(predictions,summaries,kernel):
    reports={};passes=[]
    for case,modes in predictions.items():
        lookup={r['hgvs_nt']:r for r in summaries[case]};reports[case]={}
        for mode,record in modes.items():
            ps=record['predictions'];variant_errors=[]
            for p in ps:
                r=lookup[p['hgvs_nt']];t=r['contrasts']['T']['center']
                variant_errors.append({'hgvs_nt':p['hgvs_nt'],'site':r['site'],'ROI':r['ROI'],'is_SNV':r['is_SNV'],'observed_T':t,
                                       'absolute_errors':{n:abs(F(t)-exact_prediction(p,n)) for n in NAMES}})
            def score(sub):
                sites=sorted({r['site'] for r in sub})
                site_errors={s:{n:sum(r['absolute_errors'][n] for r in sub if r['site']==s)/sum(r['site']==s for r in sub) for n in NAMES} for s in sites}
                return {'variants':len(sub),'sites':len(sites),'site_MAE':{n:sum(v[n] for v in site_errors.values())/len(sites) if sites else None for n in NAMES},'per_site_errors':site_errors}
            scores=score(variant_errors);sensitivity={'SNV':score([r for r in variant_errors if r['is_SNV']]),
                 'ROI_1':score([r for r in variant_errors if r['ROI']==1]),'ROI_2':score([r for r in variant_errors if r['ROI']==2])}
            if scores['sites']:
                criterion=kernel.expansion_criterion(scores['site_MAE']['RNA_load_linear'],[scores['site_MAE'][n] for n in NAMES[:-1]],scores['variants'],scores['sites'])
            else:criterion={'pass':False,'coverage_pass':False,'relative_comparison_defined':False,'reason':'no_scorable_predictions'}
            signs={f'{a}/{b}':0 for a in ('-1','0','1','unknown') for b in ('-1','0','1','unknown')}
            for r in summaries[case]:
                a=r['contrasts']['R']['resolved_sign'];b=r['contrasts']['T']['resolved_sign']
                signs[f'{a if a is not None else "unknown"}/{b if b is not None else "unknown"}']+=1
            reports[case][mode]={'scores':scores,'sensitivities':sensitivity,'variant_errors':variant_errors,
                 'criterion':criterion,'unavailable':record['unavailable'],'RNA_T_sign_confusion':signs,
                 'fixed_denominators':{'all_variants':125,'all_sites':53,'SNVs':113,'SNV_sites':52},
                 'mid_is_not_a_primary_rescue':mode=='mid'}
            if mode=='heavy':passes.append(criterion['pass'])
    return serial({'schema':'COMT_proxy_evaluation_v1','cases':reports,'primary_operational_expansion_pass':len(passes)==2 and all(passes),
        'primary_rule':'Coverage and >=10% reduction against ALL FOUR baselines under EACH grouping map',
        'units':'nats','inference':'Named predictive algorithms for dependent measured proxies; no calibrated uncertainty, assay equivalence or mechanism claim'})
