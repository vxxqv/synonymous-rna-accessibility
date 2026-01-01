"""Frozen full-domain FAS count analysis. No acquisition occurs in this module."""
from fractions import Fraction as F
from pathlib import Path
from statistics import median
import hashlib
import itertools
import json
import re
from ordinal_geometry import features, order_certificate, bounded_pairs

WT='GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG'
ANCESTOR='GATCCAGATCTAACTTGCTGTGGTTGTGTCTCCTGCTTCTCCCGATTCTAGTAATTGTTTGGG'
PATTERN=re.compile('GATCCAGATCTAACTTG[CGT][GT]GTGG[CT]T[GT]TGTCT[CT]CT[GT]CTT[CT]T[GC]CC[AG]ATTC[CT]A[GC]TAATTGTTTGGG')

def parse_counts(raw, expected_git_blob=None):
    if expected_git_blob:
        git_sha=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        if git_sha != expected_git_blob: raise ValueError('Source Git blob mismatch')
    lines=raw.decode('utf-8-sig').splitlines()
    if not lines or lines[0].split()!=['Genotype','Occurrences','MutationCount']:
        raise ValueError('Unexpected count header')
    out={}
    for line in lines[1:]:
        if not line.strip(): continue
        fields=line.split()
        if len(fields)!=3: raise ValueError('Malformed count row')
        s,c,m=fields
        if not PATTERN.fullmatch(s) or s in out: raise ValueError('Unexpected or duplicate genotype')
        if not c.isascii() or not c.isdigit() or not m.isascii() or not m.isdigit():
            raise ValueError('Counts/mutation number must be nonnegative integers')
        if int(m)!=sum(a!=b for a,b in zip(s,ANCESTOR)):
            raise ValueError('Ancestor mutation annotation mismatch')
        out[s]=int(c)
    if not out or sum(out.values())<=0: raise ValueError('Empty count library')
    return out

def ratio(out,inp,s):
    required=(out.get(s),out.get(WT),inp.get(s),inp.get(WT))
    if any(v is None for v in required): return None
    a,b,c,d=required
    if b<=0 or c<=0 or d<=0:return None
    return F(a*d,b*c)

def source_relative(outputs,inputs,sequences):
    """Match author's TR1-aligned normalized frequency/median-input ES then divide WT."""
    universe=set(inputs[0])
    all_libs=inputs+outputs
    aligned=[sum(v for s,v in lib.items() if s in universe) for lib in all_libs]
    if any(t<=0 for t in aligned): raise ValueError('Empty author-aligned source universe')
    def med_input(s):
        vals=[F(lib[s],aligned[k]) for k,lib in enumerate(inputs) if s in lib]
        return median(vals) if vals else None
    med={s:med_input(s) for s in set(sequences)|{WT}}
    result={s:[] for s in sequences}
    for k,lib in enumerate(outputs):
        for s in sequences:
            if s not in universe or WT not in universe or s not in lib or WT not in lib or lib[WT]<=0 or med[s] is None or med[WT] is None or med[s]<=0 or med[WT]<=0:
                result[s].append(None)
            else:
                result[s].append(F(lib[s],lib[WT])*med[WT]/med[s])
    return result,{'TR1_aligned_totals':aligned,'whole_file_totals':[sum(z.values()) for z in all_libs],
                   'total_order':'TR1..TR3, BR1..BR9','aligned_genotype_count':len(universe)}

def encode(v): return None if v is None else str(v)

def analyze(outputs,inputs,sequences):
    if len(outputs)!=9 or len(inputs)!=3:raise ValueError('Requires nine outputs and three technical inputs')
    if len(sequences)!=72 or len(set(sequences))!=72 or WT not in sequences:raise ValueError('Fixed full domain required')
    codon_starts=(23,38,29,32,41)
    states=[tuple(s[k:k+3] for k in codon_starts) for s in sequences]
    x=features(states)
    nucleotide_states=[tuple(s[k] for k in range(63) if len({q[k] for q in sequences})>1) for s in sequences]
    xn=features(nucleotide_states)
    edges=[(i,j) for i,j in itertools.combinations(range(72),2) if sum(a!=b for a,b in zip(sequences[i],sequences[j]))==1]
    if len(edges)!=204 or len(x[0])!=8:raise ValueError('Graph/features differ from verified contract')
    scenarios=[[ratio(out,inp,s) for s in sequences] for out in outputs for inp in inputs]
    lower=[];upper=[]
    for k in range(72):
        vals=[row[k] for row in scenarios]
        if any(v is None for v in vals):lower.append(F(0));upper.append(None)
        else:lower.append(min(vals));upper.append(max(vals))
    primary=order_certificate(x,lower,upper)
    nuc=order_certificate(xn,lower,upper)
    bounded=bounded_pairs(lower,upper,edges)
    source,norm=source_relative(outputs,inputs,sequences)
    point=[sum(source[s])/9 if all(v is not None for v in source[s]) else None for s in sequences]
    # Check enclosure only for complete, positive-input primary scenarios.
    for k,v in enumerate(point):
        if v is not None and upper[k] is not None and not lower[k]<=v<=upper[k]:
            raise ArithmeticError('Author-relative point lies outside scenario envelope')
    secondary=[]
    scenario_pairs={}
    for r,row in enumerate(scenarios+[point]):
        l=[F(0) if v is None else v for v in row]
        u=row
        cert=order_certificate(x,l,u)
        pair=bounded_pairs(l,u,edges)
        label=f'BR{r//3+1}_TR{r%3+1}' if r<27 else 'mean_author_per_output_relative_ES'
        scenario_pairs[label]=pair
        secondary.append({'scenario':label,'unknown_vertices':sum(v is None for v in row),
                          'codon_null':cert,'bounded_pair_count':sum(z['certified_nonempty'] for z in pair)})
    row_records=[]
    for k,s in enumerate(sequences):
        row_records.append({'vertex':k,'sequence':s,'codon_states':states[k],
                            'output_counts':[out.get(s) for out in outputs],
                            'input_counts':[inp.get(s) for inp in inputs],
                            'scenario_ratios':[encode(row[k]) for row in scenarios],
                            'lower':encode(lower[k]),'upper':encode(upper[k]),
                            'author_relative_by_output':[encode(v) for v in source[s]],
                            'point':encode(point[k])})
    return {'status':'completed','scope':'Descriptive relative included-RNA enrichment; no calibrated confidence intervals or independent validation',
            'vertices':72,'edges':204,'finite_envelopes':sum(v is not None for v in upper),
            'unknown_envelopes':sum(v is None for v in upper),'codon_design_columns':len(x[0]),
            'nucleotide_design_columns':len(xn[0]),'primary_codon_null':primary,
            'nucleotide_null_sensitivity':nuc,'conservative_barrier_pairs':sum(z['certified_nonempty'] for z in bounded),
            'normalization':norm,'vertices_data':row_records,'bounded_pairs':bounded,
            'scenario_results':secondary,'all_scenario_pairs':scenario_pairs,
            'limitations':['Same-library transfection variation only.','Empirical envelopes are not confidence intervals.',
                           'Unknown genotype-dependent recovery can mimic interactions.',
                           'Domain excludes unmeasured synonymous bypasses.',
                           'Null rejection and barriers do not attribute a particular barrier to interactions.']}

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--directory',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    freeze=json.loads((here/'freeze.json').read_text())
    for name,sha in freeze['frozen_files'].items():
        if hashlib.sha256((here/name).read_bytes()).hexdigest()!=sha:raise RuntimeError(f'Frozen file changed: {name}')
    cfg=json.loads((here/'source_contract.json').read_text())
    libs={}; receipts=[]
    for item in cfg['files']:
        p=args.directory/item['local_name']; raw=p.read_bytes()
        libs[item['label']]=parse_counts(raw,item['git_blob'])
        receipts.append({'label':item['label'],'git_blob':item['git_blob'],'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'observed_rows':len(libs[item['label']])})
    result=analyze([libs[f'BR{i}'] for i in range(1,10)],[libs[f'TR{i}'] for i in range(1,4)],cfg['sequences'])
    result['source_receipts']=receipts
    result['freeze_sha256']=hashlib.sha256((here/'freeze.json').read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','vertices','edges','finite_envelopes','unknown_envelopes','conservative_barrier_pairs')}|{'primary_status':result['primary_codon_null']['status'],'output':str(args.output)}))

if __name__=='__main__':main()
