"""Read-only alternate verification of the metadata split and coordinate labels."""
import collections
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from Bio.Seq import Seq

HERE=Path(__file__).resolve().parent
INPUT=HERE.parent/'cambray_eligibility'

def main():
    manifest=json.loads((HERE/'split_manifest.json').read_text())
    rows={r['ID']:r for r in csv.DictReader((INPUT/'sequence_metadata.tsv').open(encoding='utf-8'),delimiter='\t')}
    selected=set()
    adjacency=collections.defaultdict(set)
    for r in csv.DictReader((INPUT/'synonymous_edges.tsv').open(encoding='utf-8'),delimiter='\t'):
        adjacency[r['left_ID']].add(r['right_ID']);adjacency[r['right_ID']].add(r['left_ID'])
    unit_ids=[];diagonals=set();all_vertices=set();peptide_arms=collections.defaultdict(set)
    for unit in manifest['units']:
        ids=unit['coordinate_order_00_01_10_11'];unit_ids.append(frozenset(ids));selected.update(ids)
        for a in range(4):
            for b in range(a+1,4):
                distance=sum(x!=y for x,y in zip(rows[ids[a]]['sequence'],rows[ids[b]]['sequence']))
                self_distance=(a^b).bit_count()
                assert distance==self_distance
                assert (ids[b] in adjacency[ids[a]]) == (self_distance==1)
        expected={frozenset((ids[0],ids[3])),frozenset((ids[1],ids[2]))}
        assert {frozenset(d) for d in unit['diagonals']}==expected
        diagonals.update(expected)
        known_order=sorted(ids)
        index=int(hashlib.sha256(('C5-corner-v1|'+'|'.join(known_order)).encode()).hexdigest(),16)%4
        assert unit['withheld_ID']==known_order[index]
        #Classification endpoints are exactly the two graph neighbors of the
        #withheld vertex inside this square, not the withheld endpoint itself.
        assert set(unit['prediction_classification_endpoints'])==set(ids)&adjacency[unit['withheld_ID']]
        seen={ids[0]};pending=[ids[0]]
        while pending:
            vertex=pending.pop()
            for neighbor in adjacency[vertex]-seen:
                seen.add(neighbor);pending.append(neighbor)
        assert len(seen)==unit['component_size']
        assert set(ids)<=seen
        all_vertices.update(seen)
        peptide=str(Seq(rows[ids[0]]['sequence']).translate())
        peptide_arms[peptide].add(unit['arm'])
        assert all(str(Seq(rows[v]['sequence']).translate())==peptide for v in seen)
        assert all(rows[v]['seed']==unit['seed'] for v in seen)
    assert len(set(unit_ids))==47 and len(diagonals)==94 and len(selected)==181 and len(all_vertices)==289
    assert all(len(arms)==1 for arms in peptide_arms.values())
    seeds=sorted({u['seed'] for u in manifest['units']},key=lambda s:hashlib.sha256(('C5-joint-v1|seed='+s).encode()).digest())
    assert len(seeds)==32
    assert manifest['seed_assignment']=={s:('development' if i<16 else 'reserved_validation') for i,s in enumerate(seeds)}
    exposed=sorted({str(Seq(r['sequence']).translate()) for r in rows.values() if r['seed']=='1'})
    examined=sorted(peptide_arms)
    #Vectorized amino-acid equality uses a different representation from the
    #author split's Python zip/count comparisons.
    encode=lambda seqs:np.array([list(s.encode('ascii')) for s in seqs],dtype=np.uint8)
    x,y=encode(examined),encode(exposed)
    matches=(x[:,None,:]==y[None,:,:]).sum(axis=2)
    assert int(matches.max())==7 and not np.any(matches>=23)
    selected_seed={str(Seq(rows[v]['sequence']).translate()):rows[v]['seed'] for v in all_vertices}
    count=0;maximum=0
    for i,peptide in enumerate(examined):
        for j in range(i+1,len(examined)):
            if selected_seed[peptide]!=selected_seed[examined[j]]:
                count+=1;maximum=max(maximum,int((x[i]==x[j]).sum()))
    assert count==851 and maximum==7
    receipt={'status':'pass','squares':47,'diagonals':94,'square_vertices':181,'component_vertices':289,
             'coordinate_and_component_checks':'Exact coordinate distances and deposited graph components independently checked',
             'relatedness_checks':{'selected_pairs':count,'seed1_pairs':int(matches.size),'maximum_matches_selected':maximum,'maximum_matches_exposed':int(matches.max())},
             'distinct_algorithm_same_analyst':True,'phenotype_values_decoded':0,
             'hashes':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ('split_manifest.json','verify_split.py')}}
    (HERE/'split_corroboration.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
