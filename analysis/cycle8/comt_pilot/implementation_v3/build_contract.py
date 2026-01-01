"""Metadata-only contract builder. Does not open the source archive."""
from pathlib import Path
import hashlib,json,re


P=Path(__file__).resolve().parent
PILOT=P.parent
GATE=PILOT.parent/'comt_gate'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def seq(a,b,step=1):return list(range(a,b+1,step))


def build():
    pins={'analysis__integrative_MAVEs.Rmd':'ac267a0a5e1f07d3eadc30a2c3643a166349f0a67259ee31c649454229a869f0',
          'references__pDEST_HC_Rec_Bxb_v2_UTR5_Flag_COMT_moxGFP_LPS.fa':'725b6d3ed5bf8a3fba0b7db5a6ee649b2138b831bed9eeab7d7989665ca4c9dd'}
    for n,h in pins.items():
        if sha(GATE/n)!=h:raise ValueError('Source identity mismatch '+n)
    source=(GATE/'analysis__integrative_MAVEs.Rmd').read_text()
    # Direct transcription of notebook lines1773-1837, checked by independent
    # source review before execution. Group identities are not assay pairing.
    roles={}
    def put(role,bio,ids):
        for i in ids:
            if i in roles:raise ValueError('Duplicate role')
            roles[i]=(role,bio)
    put('NC',None,[40,108,109])
    for bio,ids in {1:[37,112,116],2:[38,39,113,117],3:[114,118],4:[115,119]}.items():put('D',bio,ids)
    for bio,ids in {1:seq(1,3)+[94,98,102],2:seq(4,6)+[95,99,103],3:[96,100,104,106],4:[97,101,105,107]}.items():put('R',bio,ids)
    for bio,ids in {1:seq(7,21)+seq(62,65)+seq(70,73)+seq(157,160)+seq(165,168),2:seq(22,36)+seq(66,69)+seq(74,77)+seq(161,164)+seq(169,172),3:seq(78,81)+seq(86,89)+seq(149,152),4:seq(82,85)+seq(90,93)+seq(153,156)}.items():put('POLY',bio,ids)
    for bio,ids in {1:seq(42,51)+[145],2:seq(52,61)+seq(140,144),3:seq(120,129)+[146],4:seq(130,139)}.items():put('FACS',bio,ids)
    rois={1:set(seq(1,46)+seq(52,56)+seq(62,69)+seq(78,85)+seq(94,97)+[104,105,108,110]+seq(112,115)+seq(120,124)+seq(130,134)+[140,146,147]+seq(149,156)+[173]),
          2:set(seq(47,51)+seq(57,61)+seq(70,77)+seq(86,93)+seq(98,103)+[106,107,109,111]+seq(116,119)+seq(125,129)+seq(135,139)+seq(141,145)+[148]+seq(157,172)+[174])}
    if rois[1]&rois[2]:raise ValueError('ROI overlap')
    bins={1:set(seq(43,58,5)+seq(121,136,5)+[141]),2:set(seq(44,59,5)+seq(122,137,5)+[140,142]),3:set(seq(45,60,5)+seq(123,143,5)),4:set(seq(46,61,5)+seq(124,144,5))}
    fractions={k:set(seq(61+k,89+k,4)+seq(148+k,168+k,4)) for k in range(1,5)}
    exclude=set([120,130,146,145,47,57,173,174]+seq(7,36))
    samples=json.loads((GATE/'sample_id_map.json').read_text())
    byid={int(r['Sample_description'].split('_')[1]):r for r in samples if re.fullmatch(r'COMT_\d+',r['Sample_description'])}
    headers=json.loads((GATE/'raw_headers.json').read_text())
    members={}; skipped=[]; train={}
    for name,header in sorted(headers.items()):
        if '/sim/' in name:
            i=1 if '108_train.' in name else 2 if '109_train.' in name else None
            if i is None:raise ValueError('Unknown simulation file')
            train[str(i)]={'summary_member':name,'truth_member':name.replace('.var.cand.vcf.summary.txt','.truth.vcf'),'header':header}
            continue
        match=re.search(r'/COMT_(\d+)_R[.]',name)
        if not match:raise ValueError('Unrecognized source member')
        i=int(match.group(1));role,bio=roles.get(i,('UNUSED',None));quadrant=next((b for b,ss in bins.items() if i in ss),None)
        frac=next((f for f,ss in fractions.items() if i in ss),None)
        use=role in ('D','R','NC') or role=='POLY' and frac in (3,4) or role=='FACS' and quadrant==3
        if i in exclude or not use:
            skipped.append({'member':name,'sample':i,'reason':'source_fixed_exclusion' if i in exclude else 'not_a_required_proxy_or_control'})
            continue
        if i not in byid:raise ValueError('Sample missing GEO metadata')
        roi=next((r for r,s in rois.items() if i in s),None)
        if roi is None:raise ValueError('Unresolved source ROI')
        title=byid[i]['Sample_title'];m=re.search(r'biorep(\d+)',title)
        geo_bio=int(m.group(1)) if m else None
        if role not in ('NC',) and geo_bio is None:raise ValueError('Missing source preparation ID')
        if role=='FACS' and not title.endswith('_P5'):raise ValueError('Low-bin title crosswalk discrepancy')
        trait='F' if role=='FACS' else ('M' if frac==3 else 'H') if role=='POLY' else role
        members[name]={'sample':'COMT_'+str(i)+'_R','sample_number':i,'ROI':roi,'trait':trait,'code_bio':bio,'geo_bio':geo_bio,'GSM':byid[i]['GSM'],'title':title,'source_quadrant':quadrant,'source_fraction':frac,'header':header}
    fasta=(GATE/list(k for k in pins if k.endswith('.fa'))[0]).read_text().splitlines()
    reference=''.join(fasta[1:]).upper();chrom=fasta[0][1:].split()[0]
    folds=json.loads((PILOT/'fold_manifest.json').read_text());catalog=json.loads((GATE/'synonymous_id_catalog.json').read_text())
    targets=[]
    for row in folds['variants']:
        edits=[(4793+p,r,a) for p,r,a in catalog['all'][row['hgvs_nt']]['edits']]
        if any(reference[p-1]!=r for p,r,a in edits):raise ValueError('Target reference mismatch')
        targets.append({**row,'edit_key':';'.join(f'{p}:{r}>{a}' for p,r,a in edits),'site':str(row['native_codon_site'])})
    required_counts={}
    for case in ('code','geo'):
        for roi in (1,2):
            for trait in ('D','R','F','H','M'):
                groups={b:sum(s['ROI']==roi and s['trait']==trait and s[case+'_bio']==b for s in members.values()) for b in range(1,5)}
                if min(groups.values())==0:raise ValueError('Missing metadata group '+str((case,roi,trait,groups)))
                required_counts[f'{case}/ROI{roi}/{trait}']=groups
    return {'schema':'COMT_A1_selective_contract_v1','status':'METADATA_ONLY_NOT_OUTCOME_AUTHORIZATION','archive_sha256':'e3fc127f268c10efbb1d1c736423c51bc5aa900e9e363d66cc8683490516cbe2','chrom':chrom,'reference':reference,'roi_positions':{'1':[4794,4898],'2':[5082,5150]},'members':members,'skipped_members':skipped,'simulation':train,'targets':targets,'expected_variants':125,'expected_sites':53,'group_sample_counts':required_counts,'input_hashes':{**pins,**{n:sha(GATE/n) for n in ('sample_id_map.json','raw_headers.json','synonymous_id_catalog.json')},'fold_manifest.json':sha(PILOT/'fold_manifest.json')},'source_executable_conventions':{'excluded_alt_codons':['TAA','TGA'],'TAG_exclusion_added':False,'source_gdna_reject_thresholds':{'1':4,'2':2},'source_caf_log10_cutoff':-5.8,'source_half_count':'.5','negative_control_absent_background':'0','low_bin':'codeQuadrant3/GEOtitleP5','background_universe':'all source-allowed variant alleles in each selected ROI position, not only synonymous targets'},'source_count_rows_read':0}


if __name__=='__main__':
    out=P/'selective_contract.json'
    if out.exists():raise FileExistsError('Fresh contract required')
    c=build();out.write_text(json.dumps(c,indent=2)+'\n')
    print(json.dumps({'members':len(c['members']),'targets':len(c['targets']),'groups':c['group_sample_counts'],'sha256':sha(out)},indent=2))
