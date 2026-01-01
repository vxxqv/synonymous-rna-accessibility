"""Guarded future integration. Reads only explicitly hash-authorized JSON/reference inputs."""
import sys
if not __debug__: raise RuntimeError('Python -O/-OO is refused before input access')
sys.dont_write_bytecode=True
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import platform
from pathlib import Path
import traceback
import warnings

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('cnot3_cohort',HERE/'cohort.py')
c=importlib.util.module_from_spec(spec); sys.modules[spec.name]=c; spec.loader.exec_module(c)
REFERENCE_PINS={'original_ids':'fac8da906ac3f2aac8c9b8ef228abdc210dc471df2cc92fb4833ce3668772e34',
                'features':'8c16aef6e62290b89d79ccf5629c61bd99d0dbb74d90e3d390046966e3d47947',
                'metadata':'5f557bbd33810668ee44114c7ace1cc0c0e9f19395f5a203846187361f931497',
                'fasta':'d29486ce127cf4ab0f681fdcfa9bcd9fe2cfd18e7331b75257eb12c4cbb385f5'}
PRINCIPAL_CHECKS={'count_5481','distinct_stable_ids','positive_finite_paired_half_lives','stored_ratio_arithmetic',
                  'positive_finite_SE_IP','STAT_IP_arithmetic','no_top1000_boundary_tie',
                  'median_log2_change_top1000_greater_rest','no_negative_SE_among_IP_candidates'}
SECONDARY_WARNINGS={'PTMD_overlap_not_available_on_full_IP_cohort','published_PTMD_cutoff_tie_rule_unreported'}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while chunk:=stream.read(1024*1024): h.update(chunk)
    return h.hexdigest()

def loads(text):
    def bad(value): raise c.IntegrityError('Nonstandard JSON numeric literal: '+value)
    def unique(pairs):
        result={}
        for key,value in pairs:
            if key in result: raise c.IntegrityError('Duplicate JSON object key: '+key)
            result[key]=value
        return result
    return json.loads(text,parse_constant=bad,object_pairs_hook=unique)

def read_json(path): return loads(Path(path).read_text(encoding='utf-8'))

def write_new(path,value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x',encoding='utf-8',newline='\n') as stream: json.dump(value,stream,indent=2,allow_nan=False); stream.write('\n')

def references(feature_path,metadata_path):
    meta={}; features={}; tids={}
    with gzip.open(metadata_path,'rt',encoding='utf-8') as stream:
        for record in map(loads,stream):
            gene=record['stable_gene_id']
            if gene in meta: raise c.IntegrityError('Duplicate selected metadata GeneID')
            meta[gene]=record
    with gzip.open(feature_path,'rt',encoding='utf-8',newline='') as stream:
        reader=csv.DictReader(stream)
        names=reader.fieldnames
        if names is None or len(names)!=len(set(names)) or set(names)!={'GeneID','transcript_id',*c.FEATURE_NAMES}: raise c.IntegrityError('Reference CSV header collision/schema')
        for row in reader:
            if None in row or any(value is None for value in row.values()): raise c.IntegrityError('Reference CSV field count')
            gene=row['GeneID']
            if gene in features: raise c.IntegrityError('Duplicate reference feature GeneID')
            if gene not in meta or meta[gene]['transcript_id']!=row['transcript_id']: raise c.IntegrityError('Reference transcript join mismatch')
            try: features[gene]={name:None if row[name]=='' else float(row[name]) for name in c.FEATURE_NAMES}
            except ValueError as exc: raise c.IntegrityError('Malformed reference numeric field') from exc
    if set(meta)!=set(features): raise c.IntegrityError('Reference metadata/features universe mismatch')
    return features,meta

def selected_proteins(fasta_path,metadata,genes,expected_sha256):
    """Select exact pinned transcript sequences, then verify all5 region hashes."""
    fasta_path=Path(fasta_path)
    if sha(fasta_path)!=expected_sha256: raise c.IntegrityError('Pinned FASTA hash mismatch')
    genes=set(genes)
    by_tid={metadata[g]['transcript_id']:g for g in genes}
    if len(by_tid)!=len(genes): raise c.IntegrityError('Transcript shared by distinct gene IDs')
    answer={}
    def check(header,sequence):
        fields=header.split('|'); tid=fields[0]
        if tid not in by_tid: return
        gene=by_tid[tid]; meta=metadata[gene]
        if gene in answer: raise c.IntegrityError('Duplicate selected FASTA transcript')
        if len(fields)<7 or c.stable_original(fields[1])!=gene or int(fields[6])!=len(sequence): raise c.IntegrityError('FASTA selected identity/length mismatch')
        regions={}
        import re
        for field in fields[7:]:
            match=re.fullmatch(r'(CDS|UTR5|UTR3):(\d+)-(\d+)',field)
            if match:
                if match[1] in regions: raise c.IntegrityError('Repeated FASTA region')
                regions[match[1]]=(int(match[2])-1,int(match[3]))
        if 'CDS' not in regions or list(regions['CDS'])!=meta['CDS_spliced_0based_halfopen']: raise c.IntegrityError('Selected CDS coordinates changed')
        begin,end=regions['CDS']
        if not 0<=begin<end<=len(sequence): raise c.IntegrityError('Invalid CDS bounds')
        if (regions.get('UTR5')!=((0,begin) if begin else None) or regions.get('UTR3')!=((end,len(sequence)) if end<len(sequence) else None)): raise c.IntegrityError('Selected UTR coordinate mismatch')
        cds=sequence[begin:end]; utr5=sequence[:begin]; utr3=sequence[end:]
        if set(sequence)-set('ACGT') or len(cds)<6 or len(cds)%3 or not cds.startswith('ATG'):
            raise c.IntegrityError('Selected complete CDS invalid')
        codons=[cds[i:i+3] for i in range(0,len(cds),3)]
        if c.CODE[codons[-1]]!='*' or any(c.CODE[x]=='*' for x in codons[:-1]): raise c.IntegrityError('Selected standard-code stop mismatch')
        protein=''.join(c.CODE[x] for x in codons[:-1])
        if len(protein)!=meta['sense_codon_count']: raise c.IntegrityError('Selected protein length mismatch')
        for name,value in [('transcript',sequence),('CDS',cds),('UTR5',utr5),('UTR3',utr3),('protein',protein)]:
            if hashlib.sha256(value.encode('ascii')).hexdigest()!=meta[name+'_sha256']: raise c.IntegrityError('Selected reference sequence hash mismatch: '+name)
        answer[gene]=protein
    with gzip.open(fasta_path,'rt',encoding='ascii') as stream:
        header=None; chunks=[]
        for line in stream:
            if line.startswith('>'):
                if header is not None: check(header,''.join(chunks))
                header=line[1:].strip(); chunks=[]
            else: chunks.append(line.strip())
        if header is not None: check(header,''.join(chunks))
    if set(answer)!=genes: raise c.IntegrityError('Selected eligible protein absent from FASTA')
    return answer

def split_payloads(eligible,features,assignments):
    """Training API gets no heldout outcome and no KO field in either predictor file."""
    by_gene={r['GeneID']:r for r in assignments}
    if len(by_gene)!=len(assignments) or set(by_gene)!={r['GeneID'] for r in eligible}: raise c.IntegrityError('Split writer ID mismatch')
    payload={'training_predictors':[],'holdout_predictors':[],'training_targets':{},'holdout_targets':{}}
    expected={'GeneID','group_ID','HL_WT','LOG2_FC_IP','SE_IP',*c.FEATURE_NAMES}
    for row in sorted(eligible,key=lambda r:r['GeneID']):
        gene=row['GeneID']; assignment=by_gene[gene]
        if assignment['outer'] not in ('training','holdout'): raise c.IntegrityError('Invalid outer split')
        predictor={'GeneID':gene,'group_ID':assignment['group_ID'],**{name:row[name] for name in ('HL_WT','LOG2_FC_IP','SE_IP')},**features[gene]}
        if set(predictor)!=expected: raise c.IntegrityError('Predictor field leak/collision')
        if not c.finite_number(row['HL_CNOT3KO']) or row['HL_CNOT3KO']<=0: raise c.IntegrityError('Nonpositive target')
        kind=assignment['outer']; payload[kind+'_predictors'].append(predictor)
        payload[kind+'_targets'][gene]=math.log(row['HL_CNOT3KO'])
    return payload

def execute(paths,authorization_path,output):
    paths={key:Path(value).resolve() for key,value in paths.items()}; output=Path(output).resolve()
    if set(paths)!={'rows','original_ids','features','metadata','fasta','upstream_gate','protocol'}: raise c.IntegrityError('Execution input roles differ')
    if output.exists() or any(path==output or output in path.parents for path in paths.values()): raise c.IntegrityError('Output must be fresh and cannot contain inputs')
    authorization=read_json(authorization_path)
    if authorization.get('operation')!='CNOT3_cohort_group_partition_v1' or authorization.get('real_execution_authorized') is not True or authorization.get('reproduction_gate_passed') is not True: raise c.IntegrityError('Explicit root execution authorization/source gate required')
    before={key:sha(path) for key,path in paths.items()}
    codes={'cohort.py':sha(HERE/'cohort.py'),'pipeline.py':sha(Path(__file__))}
    if authorization.get('input_sha256')!=before or authorization.get('code_sha256')!=codes: raise c.IntegrityError('Authorization input/code hashes differ')
    if any(before[key]!=value for key,value in REFERENCE_PINS.items()): raise c.IntegrityError('Frozen reference/source inventory changed')
    gate=read_json(paths['upstream_gate'])
    if gate.get('schema')!='CNOT3_principal_reproduction_gate_v1' or gate.get('status')!='passed' or gate.get('cohort_gene_count')!=5481 or gate.get('cohort_rows_sha256')!=before['rows']: raise c.IntegrityError('Hash-pinned upstream reproduction-gate adapter failed')
    if not isinstance(gate.get('principal_checks'),dict) or set(gate['principal_checks'])!=PRINCIPAL_CHECKS or any(value is not True for value in gate['principal_checks'].values()): raise c.IntegrityError('Principal reproduction check missing/failed')
    if gate.get('principal_discrepancies',[])!=[] or not isinstance(gate.get('secondary_warnings',[]),list) or any(w not in SECONDARY_WARNINGS for w in gate.get('secondary_warnings',[])): raise c.IntegrityError('Unknown/principal upstream discrepancy')
    if gate.get('prediction_protocol_sha256')!=before['protocol']: raise c.IntegrityError('Upstream gate refers to another prediction protocol')
    rows=read_json(paths['rows'])
    if not isinstance(rows,list): raise c.IntegrityError('Expected normalized reproduction row array')
    if any(not isinstance(r,dict) or set(r)!=c.ROW_FIELDS or not isinstance(r['GeneID'],str) for r in rows): raise c.IntegrityError('Malformed reproduction row')
    if any(not c.finite_number(r['SE_IP']) or r['SE_IP']<=0 or not c.finite_number(r['LOG2_FC_IP']) for r in rows): raise c.IntegrityError('Reproduction rows contradict positive-SE/finite-IP source gate')
    ids=sorted(r['GeneID'] for r in rows)
    if len(ids)!=5481 or len(set(ids))!=5481 or gate.get('cohort_GeneIDs')!=ids: raise c.IntegrityError('Upstream principal GeneID universe mismatch')
    import Bio
    import numpy
    environment={'implementation':platform.python_implementation(),'python':platform.python_version(),'numpy':numpy.__version__,'Biopython':Bio.__version__}
    if environment!={'implementation':'CPython','python':'3.12.14','numpy':'2.5.4','Biopython':'1.88'}: raise c.IntegrityError('Frozen initial runtime differs')
    output.mkdir(parents=True)
    write_new(output/'STARTED.json',{'status':'authorized_started','input_sha256_before':before,'code_sha256':codes,'authorization_sha256':sha(authorization_path),'environment':environment})
    def checkpoint(value):
        temporary=output/'grouping_progress.tmp'
        temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        temporary.replace(output/'grouping_progress.json')
    try:
        source=read_json(paths['original_ids'])
        if not isinstance(source,dict) or not isinstance(source.get('entries'),list): raise c.IntegrityError('Expected whitelisted original GeneID-only inventory')
        original_ids=[r['id'] for r in source['entries']]
        features,meta=references(paths['features'],paths['metadata'])
        eligible,elig_report=c.eligibility(rows,original_ids,features,meta)
        write_new(output/'eligibility.json',elig_report)
        if len(rows)!=5481 or not elig_report['reference_coverage_pass']:
            receipt={'status':'insufficient_design_coverage','sourcegate_passed':True,'input_count_5481':len(rows)==5481,
                     'reference_coverage_pass':elig_report['reference_coverage_pass'],'grouping_run':False,
                     'model_inputs_written':False,'input_sha256_before':before,'code_sha256':codes}
        else:
            proteins=selected_proteins(paths['fasta'],meta,[r['GeneID'] for r in eligible],REFERENCE_PINS['fasta'])
            with warnings.catch_warnings():
                warnings.simplefilter('error')
                grouping=c.group_proteins(proteins,progress=checkpoint)
            write_new(output/'grouping.json',grouping)
            assignments=c.assign_partitions(grouping['components'])
            write_new(output/'assignments.json',assignments)
            guards=c.design_guards(elig_report,assignments)
            write_new(output/'design_guards.json',guards)
            if guards['status']=='design_coverage_passed':
                payload=split_payloads(eligible,features,assignments)
                write_new(output/'model_inputs/training_predictors.json',payload['training_predictors'])
                write_new(output/'model_inputs/holdout_predictors.json',payload['holdout_predictors'])
                write_new(output/'model_inputs/training_targets.json',payload['training_targets'])
                write_new(output/'sealed_evaluation/holdout_targets.json',payload['holdout_targets'])
            receipt={'status':'cohort_group_partition_ready' if guards['status']=='design_coverage_passed' else guards['status'],
                     'sourcegate_passed':True,'final_eligible_genes':len(eligible),
                     'coverage_denominator':elig_report['original_membership_denominator'],
                     'grouping_complete':True,'model_inputs_written':guards['status']=='design_coverage_passed',
                     'input_sha256_before':before,'code_sha256':codes,
                     'training_API_receives_heldout_targets':False,'upstream_reproduction_saw_all_cohort_outcomes':True}
        after={key:sha(path) for key,path in paths.items()}
        if after!=before: raise c.IntegrityError('Inputs changed')
        receipt.update({'input_sha256_after':after,'inputs_immutable':True,'output_sha256':{p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*')) if p.is_file()},'fitting_run':False})
        write_new(output/'receipt.json',receipt)
        return receipt
    except BaseException as exc:
        failure={'status':'resource_inconclusive' if isinstance(exc,c.ResourceLimit) else 'integrity_failed',
                 'exception':repr(exc),'traceback':traceback.format_exc(),'input_sha256_before':before,
                 'input_sha256_after':{key:sha(path) for key,path in paths.items()},'code_sha256':codes,
                 'model_input_files_present':[p.relative_to(output).as_posix() for p in sorted(output.rglob('*')) if p.is_file() and ('model_inputs' in p.parts or 'sealed_evaluation' in p.parts)],
                 'model_inputs_usable':False,'partial_components_usable':False}
        if isinstance(exc,c.ResourceLimit): failure['grouping_receipt']=exc.receipt
        write_new(output/'FAILURE.json',failure)
        raise

def main():
    p=argparse.ArgumentParser()
    for name in ('rows','original_ids','features','metadata','fasta','upstream_gate','protocol','authorization','out'): p.add_argument('--'+name.replace('_','-'),required=True)
    args=vars(p.parse_args()); auth=args.pop('authorization'); out=args.pop('out')
    receipt=execute(args,auth,out)
    print(json.dumps({'status':receipt['status']}))
    if receipt['status']!='cohort_group_partition_ready': raise SystemExit(2)

if __name__=='__main__': main()
