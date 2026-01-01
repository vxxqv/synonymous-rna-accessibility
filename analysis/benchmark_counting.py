"""Bounded exact-arithmetic and compiled-enumerator verification cases."""
from pathlib import Path
import sys,json,time,random,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from accessibility import *
from factor_count import count_compatible
rows=[]
for m in (1,2,5,10,20,40):
    peptide='Q'*m+'L'*m;pairs=register_pairs(3*m)
    for paired in (False,True):
        t=time.perf_counter();r=count_compatible(peptide,pairs,paired_copy=paired,table_cap=200000)
        assert int(r['count'])==3**(m*(2 if paired else 1))
        rows.append({'case':'QQ_LL_aligned','codons_per_arm':m,'paired_copy':paired,'seconds':time.perf_counter()-t,**r})
rng=random.Random(1729)
for codons in (12,24,48):
    peptide=''.join(rng.choice('LRGASP') for _ in range(codons));seq=''.join(rng.choice(CODONS[a]) for a in peptide)
    remaining=list(range(len(seq)));rng.shuffle(remaining);pairs=[]
    while remaining:
        i=remaining.pop();choices=[j for j in remaining if seq[i]+seq[j] in PAIRS]
        if choices:
            j=rng.choice(choices);remaining.remove(j);pairs.append((min(i,j),max(i,j)))
    pairs=sorted(pairs)
    for paired in (False,True):
        t=time.perf_counter();r=count_compatible(peptide,pairs,paired_copy=paired,table_cap=200000)
        rows.append({'case':'seeded_matching','codons':codons,'paired_copy':paired,'seconds':time.perf_counter()-t,**r})
(ROOT/'results/counting_benchmark.json').write_text(json.dumps(rows,indent=2))
binary=ROOT/'build/exact_kernel.exe'
if not binary.exists():binary=ROOT/'build/exact_kernel'
if binary.exists():
    for _ in range(100):
        peptide=''.join(rng.choices(list(CODONS),k=2));sites=list(range(6));rng.shuffle(sites);pairs=[tuple(sorted((sites.pop(),sites.pop()))) for _ in range(rng.randrange(4))];pairs=sorted(pairs)
        process=subprocess.run([str(binary),peptide,*[f'{i}:{j}' for i,j in pairs]],check=True,capture_output=True,text=True)
        observed={line.split('\t')[0] for line in process.stdout.splitlines()}
        assert observed==set(enumerate_space(peptide,pairs))
    (ROOT/'results/compiled_kernel_validation.json').write_text(json.dumps({'cases':100,'seed':1729,'status':'all exact genotype sets agree'}))
print('exact counting and resource gates verified')
