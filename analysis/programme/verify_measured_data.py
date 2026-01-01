"""Independent finite CDS/hash replay and outcome-exposure/eligibility checks."""
import argparse
import csv
import hashlib
import gzip
import itertools
import json
from pathlib import Path
import re
from measured_data_audit import xlsx_rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    repo=Path(__file__).resolve().parents[2];audit=repo/'results/programme/measured_data_audit_v1'
    rows=list(csv.DictReader((audit/'synonymous_group_mapping.csv').open(encoding='utf-8',newline='')))
    code=dict(zip((''.join(t) for t in itertools.product('TCAG',repeat=3)), 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
    assert len(code)==64 and sum(v!='*' for v in code.values())==61
    refs={}
    for g in ['GFP','URA3']:
        text=(a.work/'primary_biology'/('Ref_'+g+'.fasta')).read_text()
        s=''.join(x.strip() for x in text.splitlines() if not x.startswith('>')).upper()
        stop=next(i+3 for i in range(0,len(s)-2,3) if code[s[i:i+3]]=='*');refs[g]=s[:stop]
    def peptide(s):return ''.join(code[s[i:i+3]] for i in range(0,len(s),3))
    def digest(s):return hashlib.sha256(s.encode()).hexdigest()
    cache={};unique=set()
    for r in rows:
        g=r['reporter'];wt=refs[g];seqs=[]
        for label in ['',r['single_A'],r['single_B'],r['group']]:
            key=(g,label)
            if key not in cache:
                s=list(wt);muts=re.findall(r'([ACGT])(\d+)([ACGT])',label)
                assert len({int(k) for _,k,_ in muts})==len(muts)
                for x,k,y in muts:
                    assert wt[int(k)-1]==x;s[int(k)-1]=y
                s=''.join(s);assert peptide(s)==peptide(wt);cache[key]=digest(s)
            seqs.append(cache[key])
        assert seqs==[r[k] for k in ['CDS_WT_sha256','CDS_A_sha256','CDS_B_sha256','CDS_AB_sha256']]
        key=(g,r['promoter'],r['medium'],r['group']);assert key not in unique;unique.add(key)
    receipt=json.loads((repo/'data/programme/measured_source_manifest.json').read_text())
    for x in receipt['files']:
        source=a.work/x['work_relative_path'];assert hashlib.sha256(source.read_bytes()).hexdigest()==x['sha256']
    supplementary=json.loads((repo/'data/programme/measured_endpoint_source_manifest.json').read_text())
    for x in supplementary['files']:
        assert hashlib.sha256((a.work/x['work_relative_path']).read_bytes()).hexdigest()==x['sha256']
    # Reconstruct source reconciliation directly, without a private work-cache derivative.
    old={d['A']['text']:d for rn,d in xlsx_rows(a.work/'primary_biology/assets/wu2022_supplement.xlsx',7)[2:] if 'A' in d}
    current=dict(xlsx_rows(a.work/'programme_sources/gfp_supplementary_data.xlsx',7))
    prior=[]
    for rn,d in current.items():
        if rn<=2:continue
        label=d['A']['text']
        if not re.fullmatch('[ACGT][0-9]+[ACGT]',label):continue
        if label in old and 'J' in d and d['J']['kind']=='n' and 'E' in old[label] and old[label]['E']['kind']=='n':
            prior.append({'genotype':label,'S7':float(d['J']['text']),'original_pool_w':float(old[label]['E']['text'])})
    assert len(prior)==2354 and sum(r['S7']==0 for r in prior)==262
    zero={r['genotype'] for r in prior if r['S7']==0}
    primary=[r for r in rows if r['column']=='J'];affected=[r for r in primary if r['single_A'] in zero or r['single_B'] in zero]
    bad=next(r for r in primary if int(r['row_AB'])==11262)
    assert bad['planned_partition_J_only']=='development' and set(bad['group'].split(','))=={'G321T','C336T'}
    assert float(current[11262]['J']['text'])==0.326463358865135
    mapping=list(csv.DictReader((repo/'data/programme/mitochondrial_sample_mapping.csv').open(encoding='utf-8',newline='')))
    geo={}
    for fname in ['mito_family.soft.gz','mito_function.soft.gz']:
        s=gzip.decompress((a.work/'primary_biology/assets'/fname).read_bytes()).decode()
        for block in s.split('^SAMPLE = ')[1:]:
            sample=block.splitlines()[0];geo[sample]=block
    assert len(mapping)==21 and all(r['sample_id'] in geo for r in mapping)
    for r in mapping:
        if r['study'].startswith('GSE173283'):
            assert '!Sample_title = '+r['label']+'\n' in geo[r['sample_id']]
    b=json.loads((audit/'structural_gate.json').read_text())
    assert len(b['primary_eligible_units'])==6 and not b['functional_outcomes_read']
    for name in ['ura3_function_v1','ura3_function_v2']:
        path=repo/'results/programme'/name;manifest=json.loads((path/'manifest.json').read_text())
        assert not (path/'data/results.json').exists() and not (path/'data/sealed_validation_predictions.json').exists()
        assert any(c['check']=='process exited successfully within enforced limits' and not c['passed'] for c in manifest['checks'])
    by_gene={}
    for g in refs:
        rr=[r for r in rows if r['reporter']==g]
        by_gene[g]={'unique_double_patterns_across_conditions':len({r['group'] for r in rr}),
                    'unique_participating_singles_across_conditions':len({r[k] for r in rr for k in ['single_A','single_B']}),
                    'unique_position_pairs_across_conditions':len({tuple(sorted([int(r['position_A']),int(r['position_B'])])) for r in rr})}
    result={'passed':True,'independent_translation_method':'Literal standard-code 64-codon table; no Biopython translation used in this replay',
            'candidate_condition_groups_verified':len(rows),'unique_reconstructed_CDS_genotypes':len(cache),'input_hashes_verified':len(receipt['files']),
            'additional_endpoint_metadata_hashes_verified':len(supplementary['files']),'sample_ids_verified':len(mapping),
            'reporters':2,'unique_patterns':by_gene,
            'primary_candidate_groups':len(primary),'known_zero_single_affected_groups':len(affected),
            'known_zero_single_affected_by_partition':{s:sum(r['planned_partition_J_only']==s for r in affected) for s in ['development','validation','quarantine']},
            'upper_bound_usable_primary_groups_given_only_known_zero_singles':len(primary)-len(affected),
            'bound_scope':'A necessary positive-finite eligibility bound only. Not the number lost to the historical <0.5 cutoff; double outcomes mostly remain unexamined.',
            'mitochondrial_primary_units':6,'mitochondrial_functional_outcomes_examined':False,
            'models_fitted':False,'validation_double_outcomes_examined':False,'frozen_failure_rules_respected':True}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
