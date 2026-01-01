"""Independent sequence/mask corroboration; input whitelist excludes phenotypes."""
import sys
if not __debug__:raise RuntimeError('Disabled assertions refused')
sys.dont_write_bytecode=True
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
from itertools import combinations

AA={'F':'UUU UUC','L':'UUA UUG CUU CUC CUA CUG','I':'AUU AUC AUA','M':'AUG',
'V':'GUU GUC GUA GUG','S':'UCU UCC UCA UCG AGU AGC','P':'CCU CCC CCA CCG',
'T':'ACU ACC ACA ACG','A':'GCU GCC GCA GCG','Y':'UAU UAC','H':'CAU CAC','Q':'CAA CAG',
'N':'AAU AAC','K':'AAA AAG','D':'GAU GAC','E':'GAA GAG','C':'UGU UGC','W':'UGG',
'R':'CGU CGC CGA CGG AGA AGG','G':'GGU GGC GGA GGG','*':'UAA UAG UGA'}
CODE={c:aa for aa,group in AA.items() for c in group.split()}
def translate(sequence):
    if len(sequence)!=96 or not set(sequence)<=set('ACGTU'):raise ValueError('Bad coding identity')
    s=sequence.replace('T','U');return ''.join(CODE[s[i:i+3]] for i in range(0,len(s),3))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preflight',required=True);parser.add_argument('--eligibility',required=True);parser.add_argument('--out',required=True);args=parser.parse_args()
    pre=Path(args.preflight);src=Path(args.eligibility);out=Path(args.out)
    if out.exists():raise FileExistsError(out)
    files=[pre/'split_manifest.json',pre/'prediction_mask_v2.json',src/'sequence_metadata.tsv',src/'components.jsonl',src/'annotated_squares.json',src/'synonymous_edges.tsv']
    before={p.name:sha(p) for p in files}
    split=json.loads(files[0].read_text());mask=json.loads(files[1].read_text())
    for name,pin in split['input_hashes'].items():
        if name in before:assert before[name]==pin
    components=[json.loads(line) for line in (src/'components.jsonl').read_text().splitlines()]
    indices={u['component_index'] for u in split['units']}
    selected={v for i in indices for v in components[i]['IDs']};rows={};exposed=set();seen=set()
    with (src/'sequence_metadata.tsv').open(encoding='utf-8',newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        assert reader.fieldnames==['ID','seed','combi','design_rep','sequence']
        for row in reader:
            assert row['ID'] not in seen;seen.add(row['ID'])
            if row['ID'] in selected:rows[row['ID']]=row
            if row['seed']=='1':exposed.add(translate(row['sequence']))
    assert len(seen)==244000 and set(rows)==selected and len(selected)==289
    peptide={v:translate(row['sequence']) for v,row in rows.items()}
    by_seed=defaultdict(set)
    for v,row in rows.items():by_seed[row['seed']].add(peptide[v])
    pairs=0;maximum=0
    seeds=sorted(by_seed,key=int)
    for i,s in enumerate(seeds):
        for t in seeds[i+1:]:
            for a in by_seed[s]:
                for b in by_seed[t]:
                    matches=sum(x==y for x,y in zip(a,b));pairs+=1;maximum=max(maximum,matches)
                    assert matches<23
    epairs=0;emax=0
    for group in by_seed.values():
        for a in group:
            for b in exposed:
                matches=sum(x==y for x,y in zip(a,b));epairs+=1;emax=max(emax,matches)
                assert matches<23
    assert len(seeds)==32 and len(exposed)==2803
    assert pairs==split['selected_cross_seed_peptide_comparisons'] and epairs==split['exposed_peptide_comparisons']
    ordered=sorted(seeds,key=lambda s:hashlib.sha256(('C5-joint-v1|seed='+s).encode()).hexdigest())
    assignment={s:('development' if i<16 else 'reserved_validation') for i,s in enumerate(ordered)}
    assert assignment==split['seed_assignment']
    edges=set()
    for a,b in combinations(sorted(selected),2):
        if peptide[a]==peptide[b] and sum(x!=y for x,y in zip(rows[a]['sequence'],rows[b]['sequence']))==1:edges.add((a,b))
    source_edges=set()
    with (src/'synonymous_edges.tsv').open(encoding='utf-8',newline='') as stream:
        for edge in csv.DictReader(stream,delimiter='\t'):
            a,b=edge['left_ID'],edge['right_ID']
            if a in selected or b in selected:
                assert a in selected and b in selected
                source_edges.add(tuple(sorted((a,b))))
    assert edges==source_edges
    for i in indices:
        nodes=set(components[i]['IDs']);reached={min(nodes)}
        while True:
            expanded=reached|{a for a,b in edges if b in reached}|{b for a,b in edges if a in reached}
            if expanded==reached:break
            reached=expanded
        assert reached==nodes
    for unit in split['units']:
        ids=unit['coordinate_order_00_01_10_11'];seq=[rows[v]['sequence'] for v in ids]
        assert len({peptide[v] for v in ids})==1 and all(assignment[rows[v]['seed']]==unit['arm'] for v in ids)
        for a,b in [(0,1),(0,2),(1,3),(2,3)]:assert tuple(sorted((ids[a],ids[b]))) in edges
        assert all(sum(x!=y for x,y in zip(seq[a],seq[b]))==2 for a,b in [(0,3),(1,2)])
        variable=[i for i,col in enumerate(zip(*seq)) if len(set(col))>1]
        assert len(variable)==2 and variable[0]//3!=variable[1]//3
        lexical=sorted(ids);target=lexical[int(hashlib.sha256(('C5-corner-v1|'+'|'.join(lexical)).encode()).hexdigest(),16)%4]
        assert target==unit['withheld_ID']
    counts={}
    for arm in ['development','reserved_validation']:
        units=[u for u in split['units'] if u['arm']==arm];union={u['withheld_ID'] for u in units}
        assert union==set(mask['arms'][arm]['masked_IDs'])
        eligible=[]
        for unit in units:
            n=len(set(unit['coordinate_order_00_01_10_11'])&union)
            candidate=next(u for u in mask['units'] if u['square_index']==unit['square_index'])
            assert candidate['prediction_eligible_from_metadata']==(n==1)
            assert candidate['mask_count']==n
            if n==1:eligible.append(unit)
        counts[arm]={'squares':len(units),'mask_IDs':len(union),'prediction_squares':len(eligible),'prediction_targets':len({u['withheld_ID'] for u in eligible}),'groups':len({u['seed'] for u in eligible})}
    assert before=={p.name:sha(p) for p in files}
    out.mkdir(parents=True)
    result={'status':'passed','input_hashes_before_after_equal':before,'records_in_sequence_only_TSV':len(seen),'selected_component_vertices':len(selected),'selected_components':len(indices),'selected_edges':len(edges),'selected_peptide_pairs':pairs,'seed1_peptide_pairs':epairs,'max_selected_matches':maximum,'max_seed1_matches':emax,'arms':counts,'phenotype_fields_or_sources_read':0,'implementation':'Handwritten standard code; direct pairwise Hamming graph; independently computed SHA split and union-mask incidence.'}
    (out/'metadata_check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
