"""Exact mitochondrial-code eligibility audit. No human variant outcomes are used."""
from pathlib import Path
from itertools import product
from collections import defaultdict
import json,argparse,sys
from Bio.Data import CodonTable
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'analysis/expansion'))
from feasibility_2sat import _codon_clauses,encode,evaluate,BITS,BASE


def run():
    data=json.loads((ROOT/'data/programme/ncbi_code2.json').read_text())
    mapping=data['RNA_mapping'];bio=CodonTable.unambiguous_dna_by_id[2]
    assert len(mapping)==64
    assert all(mapping[c.replace('T','U')]==aa for c,aa in bio.forward_table.items())
    assert all(mapping[c.replace('T','U')]=='*' for c in bio.stop_codons)
    relations=defaultdict(list)
    for codon,aa in mapping.items():
        if aa!='*':relations[aa].append(codon)
    rows=[];checks=majority_checks=0
    for aa,codons in sorted(relations.items()):
        clauses=_codon_clauses(codons);truth=set(codons)
        for x in product('ACGU',repeat=3):
            seq=''.join(x);assert evaluate(clauses,encode(seq))==(seq in truth);checks+=1
        enc={encode(c) for c in codons}
        for a,b,c in product(enc,repeat=3):
            majority=tuple(int(x+y+z>=2) for x,y,z in zip(a,b,c))
            assert majority in enc;majority_checks+=1
        rows.append({'amino_acid':aa,'codons':sorted(codons),'clauses':clauses,'clause_count':len(clauses)})
    standard={c.replace('T','U'):a for c,a in CodonTable.unambiguous_dna_by_id[1].forward_table.items()}
    changes=[{'codon':c,'standard':standard.get(c,'*'),'code2':aa} for c,aa in sorted(mapping.items()) if aa!=standard.get(c,'*')]
    assert len(relations)==20 and sum(map(len,relations.values()))==60
    assert {x['codon'] for x in changes}=={'AGA','AGG','AUA','UGA'}
    assert all(mapping[c]=='K' for c in ('AAA','AAG')) and all(mapping[c]=='F' for c in ('UUU','UUC'))
    return {'scope':'Exact finite proof of specified codon relations and compiler eligibility, not a universal reconfiguration proof',
            'table':2,'sense_codons':60,'amino_acid_relations':20,'truth_assignments_checked':checks,
            'majority_triples_checked':majority_checks,'changes_from_standard':changes,
            'maximum_clause_count':max(x['clause_count'] for x in rows),'relations':rows,
            'survives':'Two-CNF feasibility encoding for complete internal code-2 codons, including identified variables and reverse-complement literal substitution. Lys/Phe embedding is unchanged for unrestricted logical pairing graphs.',
            'requires_proof_transfer':'Component/geodesic/event-order statements require the explicit local-factor proof, not only Boolean satisfiability. Existing standard-code numerical register rates and cohort numbers are not transferred.',
            'not_modeled':['Initiator reassignment','Polyadenylation-completed stop codons','Processing boundaries','Terminal recoding','Functional effects of unmeasured variants','General reachability in overlapping frames']}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True);path=Path(parser.parse_args().out)
    result=run();path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('relations','not_modeled')}))
