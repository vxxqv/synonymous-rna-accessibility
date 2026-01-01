"""All frozen retained sample/position depth conflicts, without phenotype scores.
Selectively decodes identities, DP and frozen Boolean row-status flags only.
CAO, CAF, q, negative-control corrections, R/T/L and model outputs are not read.
"""
from pathlib import Path
from collections import defaultdict,Counter
from fractions import Fraction
import csv,gzip,hashlib,json

P=Path(__file__).resolve().parent;BASE=P.parent/'comt_results/actual_v1'
C=P.parent/'comt_pilot/implementation_v2/selective_contract.json'
FIELDS=('#CHROM','POS','REF','ALT','POS_NT','REF_NT','ALT_NT','DP','Sample','ROI','record_id')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def median(values):
    a=sorted(values);n=len(a)
    return Fraction(a[n//2]) if n%2 else Fraction(a[n//2-1]+a[n//2],2)
def selected_tsv(path,fields):
    with path.open('rb') as f:
        header=f.readline().rstrip(b'\r\n').decode('utf8').split('\t');ix={x:header.index(x) for x in fields}
        for number,line in enumerate(f,2):
            b=line.rstrip(b'\r\n').split(b'\t')
            if len(b)!=len(header):raise ValueError('Row shape')
            yield {x:b[j].decode('utf8') for x,j in ix.items()}
def edit_key(rows):
    r=rows[0];pos=int(r['POS'])
    if ',' not in r['REF']:
        edits=[(pos+i,a,b) for i,(a,b) in enumerate(zip(r['REF'],r['ALT'])) if a!=b]
    else:edits=sorted({(int(x['POS_NT']),x['REF_NT'],x['ALT_NT']) for x in rows if x['REF_NT']!=x['ALT_NT']})
    return ';'.join(f'{p}:{a}>{b}' for p,a,b in edits)

def main():
    if (P/'receipt.json').exists():raise FileExistsError('Fresh diagnosis required')
    contract=json.loads(C.read_text(encoding='utf8'));targets={r['edit_key']:r for r in contract['targets']}
    samples={r['sample']:r for r in contract['members'].values()}
    manifests={name:json.loads((BASE/name/'sealed_manifest.json').read_text(encoding='utf8')) for name in ('extract','rf')}
    if sha(C)!=manifests['extract']['input_hashes']['contract']:raise ValueError('Contract pin')
    if manifests['rf']['input_hashes']['extract_manifest']!=sha(BASE/'extract/sealed_manifest.json'):raise ValueError('Stage linkage')
    inputs={};grouped=defaultdict(list);seen=set();retained=0;rows_total=0
    for roi in (1,2):
        rp=BASE/'extract'/f'rows_{roi}.tsv';sp=BASE/'rf'/f'row_status_{roi}.tsv'
        for phase,file in [('extract',rp),('rf',sp)]:
            h=sha(file)
            if h!=manifests[phase]['files'][file.name]:raise ValueError('Input drift')
            inputs[str(file.relative_to(BASE))]=h
        statuses={}
        for s in selected_tsv(sp,('record_id','VAR_ID','frequency_pass','globally_rejected')):
            if s['record_id'] in statuses:raise ValueError('Duplicate status')
            if {s['frequency_pass'],s['globally_rejected']}-{'TRUE','FALSE'}:raise ValueError('Boolean encoding')
            statuses[s['record_id']]=s
        for row in selected_tsv(rp,FIELDS):
            rows_total+=1;rid=row['record_id']
            if rid in seen:raise ValueError('Duplicate extracted identity')
            seen.add(rid);s=statuses.pop(rid)
            ident=':'.join(row[x] for x in ('#CHROM','POS','REF','ALT'))
            if s['VAR_ID']!=ident:raise ValueError('Status identity mismatch')
            if s['frequency_pass']!='TRUE' or s['globally_rejected']=='TRUE':continue
            if not row['DP'].isdigit():raise ValueError('Depth domain')
            retained+=1;grouped[(row['Sample'],ident)].append(row)
        if statuses:raise ValueError('Unmatched status')
    positions=defaultdict(list);samebase=defaultdict(set)
    for (sample,ident),rows in sorted(grouped.items()):
        if len({r['POS'] for r in rows})!=1:raise ValueError('Variant start mismatch')
        dp=median([int(r['DP']) for r in rows]);key=edit_key(rows);target=targets.get(key)
        for r in rows:samebase[(sample,int(r['POS_NT']))].add(int(r['DP']))
        v={'VAR_ID':ident,'median_DP':str(dp),'POS':int(rows[0]['POS']),
           'raw_ref':rows[0]['REF'],'raw_alt':rows[0]['ALT'],'per_base_records':len(rows),
           'per_base_positions':[int(r['POS_NT']) for r in rows],
           'per_base_depths':[int(r['DP']) for r in rows],
           'source_record_ids':[r['record_id'] for r in rows],
           'target_hgvs':target['hgvs_nt'] if target else None,
           'single_nucleotide_raw_allele':len(rows[0]['REF'])==len(rows[0]['ALT'])==1}
        positions[(sample,v['POS'])].append(v)
    conflicts=[];direct=[];summaries=[];cause=Counter();bytrait=Counter();byroi=Counter()
    for (sample,pos),vs in sorted(positions.items()):
        depths=sorted({Fraction(v['median_DP']) for v in vs})
        if len(depths)<2:continue
        meta=samples[sample];snv=[v for v in vs if v['single_nucleotide_raw_allele']];mnp=[v for v in vs if not v['single_nucleotide_raw_allele']]
        intrinsic=any(len(set(v['per_base_depths']))>1 for v in vs)
        same_position_discordance=any(len(samebase[(sample,p)])>1 for v in vs for p in v['per_base_positions'])
        labels=[]
        if mnp:labels.append('includes_multinucleotide_raw_alleles')
        if intrinsic:labels.append('within_allele_per_base_DP_variation')
        if len({Fraction(v['median_DP']) for v in snv})>1:labels.append('SNV_depth_disagreement_at_common_start')
        if same_position_discordance:labels.append('DP_disagreement_at_same_per_base_position')
        if not labels:labels.append('distinct_medians_without_detected_within_allele_variation')
        for label in labels:cause[label]+=1
        bytrait[meta['trait']]+=1;byroi[str(meta['ROI'])]+=1
        target_ids=sorted({v['target_hgvs'] for v in vs if v['target_hgvs']})
        record={'Sample':sample,'ROI':meta['ROI'],'trait':meta['trait'],'POS':pos,
          'distinct_median_DP':[str(x) for x in depths],'depth_spread':str(depths[-1]-depths[0]),
          'collapsed_alleles':len(vs),'raw_records':sum(v['per_base_records'] for v in vs),'SNV_alleles':len(snv),'multi_alleles':len(mnp),
          'target_ids':target_ids,'causal_description_flags':labels,'same_per_base_position_DP_discordance':same_position_discordance,
          'source_order_first_VAR_ID':vs[0]['VAR_ID'],'source_order_first_median_DP':vs[0]['median_DP'],'variants':vs}
        conflicts.append(record)
        summaries.append({k:json.dumps(v) if isinstance(v,list) else v for k,v in record.items() if k!='variants'})
        for target in target_ids:direct.append({'Sample':sample,'ROI':meta['ROI'],'trait':meta['trait'],'POS':pos,'hgvs_nt':target,
               'effect_scope':'NC_ratio_input_for_matching_target' if meta['trait']=='NC' else 'direct_target_sample_denominator'})
    with (P/'all_conflicts.csv').open('x',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]) if summaries else ['Sample','POS']);w.writeheader();w.writerows(summaries)
    with gzip.open(P/'all_conflict_variant_lineage.json.gz','xt',encoding='utf8') as f:json.dump(conflicts,f,indent=2)
    with (P/'affected_target_sample_rows.csv').open('x',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=['Sample','ROI','trait','POS','hgvs_nt','effect_scope']);w.writeheader();w.writerows(direct)
    targets_direct=sorted({r['hgvs_nt'] for r in direct if r['trait']!='NC'});targets_nc=sorted({r['hgvs_nt'] for r in direct if r['trait']=='NC'})
    report={'status':'SOURCE_CONTRACT_DIAGNOSIS_ONLY_NO_CONTRASTS','raw_rows_inspected_for_identity_and_DP':rows_total,
      'rows_retained_by_unchanged_frozen_flags':retained,'collapsed_alleles':len(grouped),'sample_position_groups':len(positions),
      'conflicting_sample_positions':len(conflicts),'samples_with_conflict':len({r['Sample'] for r in conflicts}),
      'physical_start_positions_with_conflict':len({r['POS'] for r in conflicts}),
      'conflicts_by_trait':dict(sorted(bytrait.items())),'conflicts_by_ROI':dict(sorted(byroi.items())),
      'cause_flag_counts_overlapping':dict(sorted(cause.items())),
      'conflicts_containing_targets':sum(bool(r['target_ids']) for r in conflicts),
      'direct_target_sample_conflict_rows':sum(r['trait']!='NC' for r in direct),'NC_target_sample_conflict_rows':sum(r['trait']=='NC' for r in direct),
      'target_ids_with_direct_sample_conflicts':targets_direct,'target_ids_with_NC_input_conflicts':targets_nc,
      'all_affected_target_ids':sorted(set(targets_direct)|set(targets_nc)),
      'unchanged_flags_recomputed':False,'actual_CAO_CAF_q_log_R_T_L_or_predictions_read_or_computed':False,
      'depth_aggregation':'exact rational median; no WT/ratio chosen for actual records',
      'input_sha256':inputs,'contract_sha256':sha(C),'code_sha256':sha(__file__),
      'outputs_sha256':{n:sha(P/n) for n in ('all_conflicts.csv','all_conflict_variant_lineage.json.gz','affected_target_sample_rows.csv')},
      'preserved_failure':{n:sha(BASE/n) for n in ('normalize_execution.log','normalize_execution_receipt.json')}}
    (P/'receipt.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('target_ids_with_direct_sample_conflicts','target_ids_with_NC_input_conflicts','all_affected_target_ids','input_sha256','outputs_sha256','preserved_failure')},indent=2))
    print('Affected target IDs',len(report['all_affected_target_ids']))

if __name__=='__main__':main()
