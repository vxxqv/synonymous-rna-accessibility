"""Adaptive positive-coverage follow-up to the mostly-UNSAT random compiler audit."""
from pathlib import Path
from itertools import product, combinations
import argparse
import json
import random

from overlap_probe import occurrences, labels, compile_overlap, direct, PAIR
from feasibility_2sat import encode, evaluate, solve_cnf, BASE


def run():
    rng = random.Random(73692)
    rows=[]
    truth_checks=0
    for k in range(20):
        n=6
        second=[('+',1),('+',2),('-',0),('-',1),('-',2)][k%5]
        codons=occurrences(n,'+',0)+occurrences(n,*second)
        for attempt in range(100):
            reference=''.join(rng.choice('ACGU') for _ in range(n))
            peptide=labels(reference,codons)
            if peptide is not None:
                break
        else:
            raise AssertionError('Positive source not found within bound')
        compatible=[e for e in combinations(range(n),2) if reference[e[0]]+reference[e[1]] in PAIR]
        edges=tuple(rng.sample(compatible,min(k%3,len(compatible))))
        fixed={k%n:reference[k%n]} if k%2 else {}
        nvars,clauses=compile_overlap(n,codons,peptide,edges,fixed)
        assert direct(reference,codons,peptide,edges,fixed)
        assert evaluate(clauses,encode(reference))
        count=0
        for seq in map(''.join,product('ACGU',repeat=n)):
            exact=direct(seq,codons,peptide,edges,fixed)
            assert exact==evaluate(clauses,encode(seq))
            truth_checks+=1
            count+=exact
        answer=solve_cnf(nvars,clauses)
        assert answer is not None and count>0
        witness=''.join(BASE[answer[2*p:2*p+2]] for p in range(n))
        assert direct(witness,codons,peptide,edges,fixed)
        rows.append({'frame':second,'reference':reference,'peptide':peptide,
                     'pair_edges':edges,'pins':fixed,'solutions':count,'solver_witness':witness})
    # Directly contradictory codon obligations and a reverse-complement contradiction.
    for complement,positions in [(False,(0,1,2)),(True,(2,1,0))]:
        nvars,clauses=compile_overlap(3,[(positions,complement),(positions,complement)],['K','F'])
        assert solve_cnf(nvars,clauses) is None
        assert not any(evaluate(clauses,encode(''.join(x))) for x in product('ACGU',repeat=3))
    return {'scope':'Adaptive implementation follow-up, not a new discovery or prospective validation',
            'known_SAT_instances':len(rows),'assignment_truth_checks':truth_checks,
            'explicit_UNSAT_instances':2,'rows':rows}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out',required=True)
    path=Path(parser.parse_args().out);result=run()
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
