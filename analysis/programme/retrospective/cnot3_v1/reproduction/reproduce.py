"""Bounded published-control reproduction. Real outcomes require a root freeze.

No model fitting. XML values are decoded only for explicitly selected columns
and rows. Imported use supports synthetic tests without opening any workbook.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import posixpath
import statistics
import sys
import xml.etree.ElementTree as ET
import zipfile

NS = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
FIELDS = {'A':'GENE_ID','B':'HL_WT','C':'HL_CNOT3KO','D':'HL_FC_LOG2',
          'F':'PTMD_SCORE_ABSOLUTE','G':'PTMD_SCORE',
          'I':'LOG2_FC_IP','J':'SE_IP','K':'STAT_IP'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def selected_cells(path, columns, rows=None, audit=None):
    """Read only TableS2 requested cells. rows=None means every existing row.

    Nonselected cell values, formulas and shared strings are not decoded. A
    selected formula is rejected even if a cached value is present.
    """
    selected = {}
    with zipfile.ZipFile(path) as z:
        book=ET.fromstring(z.read('xl/workbook.xml'))
        rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
        mapping={x.attrib['Id']:x.attrib['Target'] for x in rels}
        sheets=[x for x in book.findall('.//s:sheet',NS) if x.attrib['name']=='TableS2']
        if len(sheets)!=1:
            raise ValueError('Expected exactly one TableS2 sheet')
        target=mapping[sheets[0].attrib['{'+NS['r']+'}id']]
        member=target.lstrip('/') if target.startswith('/') else posixpath.normpath(posixpath.join('xl',target))
        if not member.startswith('xl/') or '..' in member.split('/'):
            raise ValueError('Invalid worksheet relationship')
        with z.open(member) as stream:
            for _, elem in ET.iterparse(stream,events=('end',)):
                if elem.tag!='{'+NS['s']+'}row':
                    continue
                row=int(elem.attrib['r'])
                if rows is None or row in rows:
                    for cell in elem.findall('s:c',NS):
                        coord=cell.attrib['r']
                        col=''.join(x for x in coord if x.isalpha())
                        if col not in columns:
                            continue
                        if cell.find('s:f',NS) is not None:
                            raise ValueError('Selected formula cell: '+coord)
                        typ=cell.attrib.get('t','n')
                        if typ=='inlineStr':
                            value=''.join(x.text or '' for x in cell.findall('.//s:t',NS))
                        else:
                            value=cell.findtext('s:v',None,NS)
                        selected[(row,col)]=(typ,value)
                elem.clear()
        needed={int(v) for t,v in selected.values() if t=='s' and v is not None}
        shared={}
        if needed:
            with z.open('xl/sharedStrings.xml') as stream:
                i=0
                for _,elem in ET.iterparse(stream,events=('end',)):
                    if elem.tag=='{'+NS['s']+'}si':
                        if i in needed:
                            shared[i]=''.join(x.text or '' for x in elem.findall('.//s:t',NS))
                        i+=1
                        elem.clear()
            if needed!=set(shared):
                raise ValueError('Missing selected shared string')
        result={}
        for (row,col),(typ,value) in selected.items():
            if typ=='s' and value is not None:
                value=shared[int(value)]
            elif typ=='e':
                value={'cell_error':value}
            result.setdefault(row,{})[FIELDS[col]]=value
        if audit is not None:
            audit.append({'columns':sorted(columns),'row_selector':'all' if rows is None else sorted(rows),
                          'decoded_cells':[str(row)+':'+col for row,col in sorted(selected)]})
        return result


def number(value):
    if isinstance(value,bool) or isinstance(value,dict):
        return None
    try:
        ans=float(value)
    except (TypeError,ValueError):
        return None
    return ans if math.isfinite(ans) else None


def ranks(values):
    order=sorted(range(len(values)),key=lambda i:values[i])
    ans=[0.0]*len(values)
    i=0
    while i<len(order):
        j=i+1
        while j<len(order) and values[order[j]]==values[order[i]]:
            j+=1
        for k in range(i,j):
            ans[order[k]]=(i+1+j)/2
        i=j
    return ans


def spearman(x,y):
    if len(x)!=len(y) or len(x)<2:
        return None
    a,b=ranks(x),ranks(y)
    ma,mb=statistics.mean(a),statistics.mean(b)
    aa=sum((v-ma)**2 for v in a); bb=sum((v-mb)**2 for v in b)
    return None if aa==0 or bb==0 else sum((u-ma)*(v-mb) for u,v in zip(a,b))/math.sqrt(aa*bb)


def ks_statistic(x,y):
    if not x or not y:
        return None
    a,b=sorted(x),sorted(y)
    return max(abs(bisect_right(a,v)/len(a)-bisect_right(b,v)/len(b)) for v in set(a+b))


def distribution(x,y,pvalue=False):
    if not x or not y:
        return {'status':'empty_group','n_high':len(x),'n_remainder':len(y)}
    b=sorted(y)
    wins=sum(bisect_left(b,v)+(bisect_right(b,v)-bisect_left(b,v))/2 for v in x)
    superiority=wins/(len(x)*len(y))
    ans={'status':'computed','n_high':len(x),'n_remainder':len(y),
         'median_high':statistics.median(x),'median_remainder':statistics.median(y),
         'median_difference':statistics.median(x)-statistics.median(y),
         'mean_difference':statistics.mean(x)-statistics.mean(y),
         'probability_high_greater_with_half_ties':superiority,
         'rank_biserial':2*superiority-1,'ks_statistic':ks_statistic(x,y),
         'ks_pvalue':None,'ks_pvalue_interpretation':'Not a biological-replicate significance test'}
    if pvalue:
        try:
            import scipy
            from scipy.stats import ks_2samp
        except ImportError:
            ans['ks_pvalue_status']='omitted_existing_runtime_unavailable'
        else:
            test=ks_2samp(x,y,alternative='two-sided',method='asymp')
            assert abs(float(test.statistic)-ans['ks_statistic'])<1e-12
            ans.update(ks_pvalue=float(test.pvalue),scipy_version=scipy.__version__,
                       ks_pvalue_status='conventional_two_sided_asymptotic_gene_distribution_reference_only',
                       ks_pvalue_caveat='Independent continuous observations are the reference assumption; ties and gene dependence can invalidate calibration. No replicate or mechanistic inference.')
    return ans


def top_group(records,field,n):
    if len(records)<n:
        return [],records,{'status':'insufficient_rows','n_available':len(records)}
    ordered=sorted(records,key=lambda r:(-r[field],r['GENE_ID']))
    cutoff=ordered[n-1][field]
    equal=[r['GENE_ID'] for r in ordered if r[field]==cutoff]
    crosses=len(ordered)>n and ordered[n][field]==cutoff
    return ordered[:n],ordered[n:],{'status':'selected','cutoff_value':cutoff,
        'boundary_tie_crosses':crosses,'all_cutoff_tie_ids':equal,
        'tie_rule':'Exact numeric equality; ascending exact GeneID for fixed membership; no tolerance tie collapse'}


def ip_candidates(rows):
    ids=Counter(r.get('GENE_ID') for r in rows if isinstance(r.get('GENE_ID'),str))
    candidates=[]; exclusions=[]; issues=[]
    for raw in rows:
        reasons=[]; gid=raw.get('GENE_ID')
        if not isinstance(gid,str) or not gid or gid.strip()!=gid:
            reasons.append('invalid_gene_id')
        elif ids[gid]!=1:
            reasons.append('duplicate_gene_id')
        record=dict(raw)
        for field in ['STAT_IP','LOG2_FC_IP','SE_IP']:
            value=number(raw.get(field))
            if value is None:
                reasons.append('missing_or_nonfinite_'+field)
            record[field]=value
        if reasons:
            exclusions.append({'row':raw.get('_row'),'GENE_ID':gid,'stage':'IP','reasons':reasons})
        else:
            candidates.append(record)
            if record['SE_IP']<0:
                issues.append('negative_SE_IP:'+gid)
    if any(c>1 for c in ids.values()):
        issues.append('duplicate_GeneIDs_not_resolved')
    if any('invalid_gene_id' in x['reasons'] for x in exclusions):
        issues.append('invalid_GeneIDs_not_resolved')
    return candidates,exclusions,issues


def analyze(rows,protocol):
    """Pure analysis for supplied records, usable with synthetic objects only."""
    candidates,excluded,issues=ip_candidates(rows)
    cohort=[]; ratio_diag=[]; stat_diag=[]
    for raw in candidates:
        row=dict(raw); reasons=[]
        for field in ['HL_WT','HL_CNOT3KO','HL_FC_LOG2']:
            row[field]=number(raw.get(field))
            if row[field] is None:
                reasons.append('missing_or_nonfinite_'+field)
        for field in ['HL_WT','HL_CNOT3KO']:
            if row[field] is not None and row[field]<=0:
                reasons.append('nonpositive_'+field)
        if reasons:
            excluded.append({'row':raw.get('_row'),'GENE_ID':raw['GENE_ID'],'stage':'HL','reasons':reasons})
            continue
        expected=math.log2(row['HL_CNOT3KO'])-math.log2(row['HL_WT'])
        residual=row['HL_FC_LOG2']-expected
        tol=protocol['ratio_atol']+protocol['ratio_rtol']*max(abs(expected),abs(row['HL_FC_LOG2']))
        ratio_diag.append({'GENE_ID':row['GENE_ID'],'computed_log2_ratio':expected,
                           'stored_minus_computed':residual,'tolerance':tol,'consistent':abs(residual)<=tol})
        if abs(residual)>tol:
            issues.append('stored_HL_ratio_inconsistent:'+row['GENE_ID'])
        if row['SE_IP']==0:
            stat_diag.append({'GENE_ID':row['GENE_ID'],'computed_stat':None,'reason':'zero_SE_IP'})
            issues.append('zero_SE_IP:'+row['GENE_ID'])
        else:
            val=row['LOG2_FC_IP']/row['SE_IP']
            if not math.isfinite(val):
                stat_diag.append({'GENE_ID':row['GENE_ID'],'computed_stat':None,'reason':'nonfinite_LFC_over_SE'})
                issues.append('nonfinite_LFC_over_SE:'+row['GENE_ID'])
            else:
                discrepancy=row['STAT_IP']-val
                stat_diag.append({'GENE_ID':row['GENE_ID'],'computed_stat':val,
                                  'archived_minus_LFC_over_SE':discrepancy,
                                  'consistent':abs(discrepancy)<=protocol['stat_atol']+protocol['stat_rtol']*max(abs(val),abs(row['STAT_IP']))})
                if not stat_diag[-1]['consistent']:
                    issues.append('archived_IP_statistic_discrepancy:'+row['GENE_ID'])
        for field in ['PTMD_SCORE','PTMD_SCORE_ABSOLUTE']:
            row[field]=number(raw.get(field))
        cohort.append(row)
    if len(cohort)!=protocol['expected_cohort_n']:
        issues.append('published_cohort_size_mismatch')
    n=protocol['top_n']
    high,rest,ties=top_group(cohort,'STAT_IP',n)
    if ties.get('boundary_tie_crosses'):
        issues.append('published_IP_cutoff_tie_rule_unreported')
    if not high or not rest:
        issues.append('insufficient_comparison_groups')
    metrics={field:distribution([r[field] for r in high if r[field] is not None],
                                [r[field] for r in rest if r[field] is not None],
                                protocol['conventional_ks_pvalue'])
             for field in ['HL_FC_LOG2','HL_WT','HL_CNOT3KO','PTMD_SCORE']}
    fields=['HL_FC_LOG2','PTMD_SCORE','PTMD_SCORE_ABSOLUTE','STAT_IP']
    complete=[r for r in cohort if all(r[f] is not None for f in fields)]
    incomplete=[{'GENE_ID':r['GENE_ID'],'missing':[f for f in fields if r[f] is None]} for r in cohort if not all(r[f] is not None for f in fields)]
    correlations={a:{b:spearman([r[a] for r in complete],[r[b] for r in complete]) for b in fields} for a in fields}
    overlap={'status':'not_computed','universe':'same reconstructed Figure2 cohort; harmonized diagnostic, not authenticated original top-PTMD universe'}
    if all(r['PTMD_SCORE'] is not None for r in cohort) and len(cohort)>=n:
        ptmd,_,pties=top_group(cohort,'PTMD_SCORE',n)
        common=sorted(set(r['GENE_ID'] for r in high)&set(r['GENE_ID'] for r in ptmd))
        overlap.update(status='computed',count=len(common),ids=common,top_ptmd_ids=[r['GENE_ID'] for r in ptmd],tie_diagnostic=pties,
                       published_count=279,chance_expected_if_same_uniform_universe=n*n/len(cohort))
        if pties.get('boundary_tie_crosses'):
            issues.append('published_PTMD_cutoff_tie_rule_unreported')
    else:
        issues.append('PTMD_overlap_not_available_on_full_IP_cohort')
    return {'status':'descriptive_reconstruction_only','n_input':len(rows),'n_IP_candidates':len(candidates),
        'n_cohort':len(cohort),'n_high':len(high),'n_remainder':len(rest),
        'expected_published_n':protocol['expected_cohort_n'],
        'cohort_ids':[r['GENE_ID'] for r in cohort],'high_ids':[r['GENE_ID'] for r in high],
        'remainder_ids':[r['GENE_ID'] for r in rest],'selection_ties':ties,'exclusions':excluded,
        'half_life_ratio_diagnostics':ratio_diag,'IP_statistic_diagnostics':stat_diag,
        'group_effects':metrics,'spearman_complete_case_n':len(complete),
        'spearman_complete_case_exclusions':incomplete,'spearman':correlations,'PTMD_overlap':overlap,
        'discrepancies':issues,'unresolved_source_rules':protocol['unresolved_source_rules'],
        'extension_gate':'HOLD until root resolves discrepancies and source-rule limits; never auto-approved',
        'interpretation':'Published-control reconstruction, not discovery or independent biological replication. Published figure axes label minutes; source-column scale and fitting calibration remain unresolved. Stored outcomes and archived STAT are retained, never replaced by diagnostics.'}


def execute(workbook,protocol_path,freeze_path,outdir):
    out=Path(outdir)
    if out.exists():
        raise FileExistsError('Output path must be new')
    freeze=json.loads(Path(freeze_path).read_text(encoding='utf-8'))
    if freeze.get('status')!='root_authorized_outcome_read':
        raise ValueError('Root outcome authorization freeze missing')
    for label,path in [('protocol',protocol_path),('implementation',__file__),('workbook',workbook)]:
        if freeze.get(label+'_sha256')!=sha(path):
            raise ValueError(label+' hash mismatch before outcome decoding')
    protocol=json.loads(Path(protocol_path).read_text(encoding='utf-8'))
    if protocol['source_sheet']!='TableS2' or protocol['selected_columns']!=FIELDS:
        raise ValueError('Unsupported frozen schema')
    if protocol['source_workbook_sha256']!=freeze['workbook_sha256']:
        raise ValueError('Protocol and freeze disagree on source workbook')
    audit=[]
    first=selected_cells(workbook,{'A','I','J','K'},audit=audit)
    for col in ['A','I','J','K']:
        if first.get(1,{}).get(FIELDS[col])!=FIELDS[col]:
            raise ValueError('Source header mismatch '+col)
    rows=[dict(data,_row=row) for row,data in sorted(first.items()) if row>1]
    candidates,_,_=ip_candidates(rows)
    candidate_rows={r['_row'] for r in candidates}
    second=selected_cells(workbook,{'B','C','D'},candidate_rows|{1},audit)
    for col in ['B','C','D']:
        if second.get(1,{}).get(FIELDS[col])!=FIELDS[col]:
            raise ValueError('Source header mismatch '+col)
    for row in rows:
        row.update(second.get(row['_row'],{}))
    paired={r['_row'] for r in rows if r['_row'] in candidate_rows and
            (number(r.get('HL_WT')) or 0)>0 and (number(r.get('HL_CNOT3KO')) or 0)>0 and number(r.get('HL_FC_LOG2')) is not None}
    third=selected_cells(workbook,{'F','G'},paired|{1},audit)
    for col in ['F','G']:
        if third.get(1,{}).get(FIELDS[col])!=FIELDS[col]:
            raise ValueError('Source header mismatch '+col)
    for row in rows:
        row.update(third.get(row['_row'],{}))
    result=analyze(rows,protocol)
    out.mkdir(parents=True)
    (out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    retained=set(result['cohort_ids'])
    normalized=[{'GeneID':r['GENE_ID'],**{k:number(r.get(k)) for k in ['HL_WT','HL_CNOT3KO','LOG2_FC_IP','SE_IP']}} for r in rows if r.get('GENE_ID') in retained]
    (out/'cohort_rows.json').write_text(json.dumps(normalized,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    (out/'exposure.json').write_text(json.dumps({'stages':audit,'no_context_features_decoded':True,
        'caveat':'All reconstructed cohort outcomes contribute to this published-control summary. Later modeling of the same endpoint is not wholly unexamined validation.'},indent=2)+'\n',encoding='utf-8')
    (out/'receipt.json').write_text(json.dumps({'input_freeze':freeze,'freeze_sha256':sha(freeze_path),
        'result_sha256':sha(out/'result.json'),'cohort_sha256':sha(out/'cohort_rows.json'),'exposure_sha256':sha(out/'exposure.json'),
        'python':sys.version,'phase':'authorized_published_reproduction','new_models_fitted':False},indent=2)+'\n',encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol',required=True)
    parser.add_argument('--workbook',required=True)
    parser.add_argument('--freeze',required=True)
    parser.add_argument('--out',required=True)
    args=parser.parse_args()
    result=execute(args.workbook,args.protocol,args.freeze,args.out)
    print(json.dumps({k:result[k] for k in ['status','n_cohort','n_high','n_remainder','extension_gate']}))


if __name__=='__main__':
    main()
