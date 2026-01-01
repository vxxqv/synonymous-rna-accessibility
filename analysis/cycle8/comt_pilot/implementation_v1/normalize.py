"""Source count normalization using exact rational count arithmetic."""
from collections import defaultdict
from fractions import Fraction as F
import math


def median(values):
    a=sorted(values)
    if not a:raise ValueError('Empty median')
    n=len(a);return a[n//2] if n%2 else (a[n//2-1]+a[n//2])/2


def raw_id(row):return ':'.join(str(row[n]) for n in ('#CHROM','POS','REF','ALT'))


def edit_key(rows):
    first=rows[0];pos=int(first['POS']);ref,alt=first['REF'],first['ALT']
    if ',' not in ref:
        edits=[(pos+i,a,b) for i,(a,b) in enumerate(zip(ref,alt)) if a!=b]
    else:
        edits=sorted({(int(r['POS_NT']),r['REF_NT'],r['ALT_NT']) for r in rows if r['REF_NT']!=r['ALT_NT']})
    return ';'.join(f'{p}:{a}>{b}' for p,a,b in edits)


def normalized(rows,statuses,contract):
    byrecord={r['record_id']:r for r in rows}
    if len(byrecord)!=len(rows) or set(byrecord)!=set(statuses):raise ValueError('Complete unique row-status alignment required')
    samples={m['sample']:m for m in contract['members'].values()}
    grouped=defaultdict(list); exclusions=[]
    for r in rows:
        s=statuses[r['record_id']]
        if s['VAR_ID']!=raw_id(r):raise ValueError('Status identity disagreement')
        if not isinstance(s['frequency_pass'],bool) or not isinstance(s['globally_rejected'],bool):raise ValueError('Boolean row-status fields required')
        if not s['frequency_pass'] or s['globally_rejected']:
            exclusions.append({'record_id':r['record_id'],'raw_id':raw_id(r),'frequency_pass':s['frequency_pass'],'globally_rejected':s['globally_rejected']})
            continue
        if r['Sample'] not in samples:raise ValueError('Unknown sample')
        grouped[(r['Sample'],raw_id(r))].append(r)
    collapsed=[]
    for (sample,ident),rr in sorted(grouped.items()):
        pos={int(r['POS']) for r in rr}
        if len(pos)!=1:raise ValueError('Variant position conflict')
        dp=median([F(str(r['DP'])) for r in rr]);cao=median([F(str(r['CAO'])) for r in rr])
        collapsed.append({'Sample':sample,'VAR_ID':ident,'POS':pos.pop(),'DP':dp,'CAO':cao,
                          'ROI':samples[sample]['ROI'],'edit_key':edit_key(rr),
                          'source_records':[r['record_id'] for r in rr]})
    positions=defaultdict(list)
    for r in collapsed:positions[(r['Sample'],r['POS'])].append(r)
    lineage=[]
    for (sample,pos),rr in sorted(positions.items()):
        dp={r['DP'] for r in rr}
        if len(dp)!=1:raise ValueError('Conflicting source-median DP at sample/position')
        depth=dp.pop();total=sum(r['CAO'] for r in rr);wt=depth-total
        if wt<0:raise ValueError('Negative reconstructed WT denominator')
        lineage.append({'Sample':sample,'POS':pos,'DP':str(depth),'sum_retained_CAO':str(total),'WT':str(wt),
                        'alleles':[{'VAR_ID':r['VAR_ID'],'CAO':str(r['CAO'])} for r in rr]})
        for r in rr:r['WT']=wt;r['q']=(r['CAO']+F(1,2))/(wt+F(1,2))
    controls=defaultdict(list)
    for r in collapsed:
        if samples[r['Sample']]['trait']=='NC':controls[(r['ROI'],r['VAR_ID'])].append(r)
    target_lookup={r['edit_key']:r for r in contract['targets']}
    if len(target_lookup)!=len(contract['targets']):raise ValueError('Duplicate target edit key')
    measurements={};target_sources=defaultdict(set)
    for r in collapsed:
        meta=samples[r['Sample']]
        if meta['trait']=='NC':continue
        nc=controls.get((r['ROI'],r['VAR_ID']),[])
        background=median([x['q'] for x in nc]) if nc else F(0)
        value=r['q']-background
        if r['edit_key'] not in target_lookup:continue
        target=target_lookup[r['edit_key']]
        if target['ROI']!=r['ROI']:raise ValueError('Target ROI mismatch')
        key=(target['hgvs_nt'],r['Sample'])
        if key in measurements:raise ValueError('Multiple raw records represent one target/sample')
        target_sources[target['hgvs_nt']].add(r['VAR_ID'])
        measurements[key]={'hgvs_nt':target['hgvs_nt'],'Sample':r['Sample'],'trait':meta['trait'],
             'q':str(r['q']),'negative_control_q':str(background),'corrected_q':str(value),
             'nc_samples':[x['Sample'] for x in nc],'nc_absence_default_used':not nc,
             'raw_id':r['VAR_ID'],'log_q':math.log(float(value)) if value>0 else None,
             'status':'finite' if value>0 else 'nonpositive_after_source_NC_subtraction'}
        if value>0 and not math.isfinite(measurements[key]['log_q']):raise ValueError('Nonfinite natural logarithm')
    all_measurements=[]
    for target in contract['targets']:
        for sample,meta in sorted(samples.items()):
            if meta['ROI']!=target['ROI'] or meta['trait']=='NC':continue
            key=(target['hgvs_nt'],sample)
            if key in measurements:all_measurements.append(measurements[key])
            else:all_measurements.append({'hgvs_nt':target['hgvs_nt'],'Sample':sample,'trait':meta['trait'],'log_q':None,'status':'absent_after_declared_source_filter_or_not_in_source'})
    return {'schema':'COMT_source_normalization_v1','units':'natural_log_corrected_variant_to_position_WT_ratio',
            'measurements':all_measurements,'denominator_lineage':lineage,'row_exclusions':exclusions,
            'target_raw_representations':{k:sorted(v) for k,v in target_sources.items()},
            'source_allowed_background_alleles':len(collapsed),'unavailable_is_not_zero':True}


def contrast_summaries(result,contract):
    samples={m['sample']:m for m in contract['members'].values()}
    lookup={(r['hgvs_nt'],r['Sample']):r for r in result['measurements']}
    out={}
    for case in ('code','geo'):
        rows=[]
        for target in contract['targets']:
            traits={}
            for trait in ('D','R','F','H','M'):
                groups=[]
                for bio in range(1,5):
                    ids=sorted(s for s,m in samples.items() if m['ROI']==target['ROI'] and m['trait']==trait and m[case+'_bio']==bio)
                    if not ids:raise ValueError('Missing structural preparation group')
                    values=[lookup[(target['hgvs_nt'],s)]['log_q'] for s in ids]
                    finite=[x for x in values if x is not None]
                    groups.append({'bio':bio,'samples':ids,'values':values,'technical_finite':len(finite),
                                   'median':float(median(finite)) if finite else None})
                finite=[g['median'] for g in groups if g['median'] is not None]
                complete=len(finite)==4
                traits[trait]={'groups':groups,'complete':complete,'center':float(median(finite)) if complete else None,
                               'lower':min(finite) if complete else None,'upper':max(finite) if complete else None}
            contrasts={}
            for name,a,b in [('R','R','D'),('T','D','F'),('L','H','R'),('L_mid','M','R')]:
                aa,bb=traits[a],traits[b];complete=aa['complete'] and bb['complete']
                value={'complete':complete,'center':aa['center']-bb['center'] if complete else None,
                       'lower':aa['lower']-bb['upper'] if complete else None,
                       'upper':aa['upper']-bb['lower'] if complete else None}
                value['resolved_sign']=None if not complete else 1 if value['lower']>0 else -1 if value['upper']<0 else 0
                contrasts[name]=value
            rows.append({**target,'assays':traits,'contrasts':contrasts,
                         'primary_complete':all(contrasts[k]['complete'] for k in ('R','T','L')),
                         'mid_complete':all(contrasts[k]['complete'] for k in ('R','T','L_mid'))})
        out[case]=rows
    return out
