"""Sequence/metadata audit and structural feature gate; never reads functional values.

External source files remain outside this repository. See measured_source_manifest.json.
"""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import tarfile
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
UNITS = {'ND1':['ND1'], 'ND2':['ND2'], 'COX1':['COX1'], 'COX2':['COX2'],
         'ATP68':['ATP8','ATP6'], 'COX3':['COX3'], 'ND3':['ND3'],
         'ND44L':['ND4L','ND4'], 'ND5':['ND5'], 'ND6':['ND6'], 'CYTB':['CYTB']}
COLS = {'B':('GFP','PTDH3','SC'), 'C':('GFP','PAGP1','SC'),
        'D':('URA3','PTDH3','SC'), 'E':('URA3','PAGP1','SC'),
        'F':('GFP','PTDH3','YPD'), 'G':('GFP','PAGP1','YPD'),
        'H':('URA3','PTDH3','YPD'), 'I':('URA3','PAGP1','YPD'),
        'J':('URA3','PTDH3','SC-URA'), 'K':('URA3','PAGP1','SC-URA')}

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n', encoding='utf-8')
def write_csv(path, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def xlsx_rows(path, sheet):
    """Return XML cells; callers choose whether numerical values may be accessed."""
    with ZipFile(path) as z:
        strings=[''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)] if 'xl/sharedStrings.xml' in z.namelist() else []
        t=ET.fromstring(z.read(f'xl/worksheets/sheet{sheet}.xml'))
        rows=[]
        for row in t.findall('s:sheetData/s:row',NS):
            d={}
            for c in row.findall('s:c',NS):
                col=re.match('[A-Z]+',c.get('r')).group();v=c.find('s:v',NS);kind=c.get('t','n')
                if v is not None:
                    d[col]={'kind':kind, 'text': strings[int(v.text)] if kind=='s' else v.text}
                elif kind=='inlineStr': d[col]={'kind':'s','text':''.join(c.find('s:is',NS).itertext())}
            rows.append((int(row.get('r')),d))
        return rows

def mutations(label):
    ms=re.findall(r'([ACGT])(\d+)([ACGT])',label)
    if re.sub(r'[ACGT]\d+[ACGT]','',label).strip(' ,;:_-+') or not ms:
        raise ValueError('Unparsed label '+label)
    if len({int(p) for a,p,b in ms})!=len(ms): raise ValueError('Duplicate site '+label)
    return tuple(sorted((a,int(p),b) for a,p,b in ms))
def mlabel(ms): return ','.join(f'{a}{p}{b}' for a,p,b in ms)

def component_sizes(pairs):
    adj=defaultdict(set)
    for a,b in pairs: adj[a].add(b);adj[b].add(a)
    done=set();sizes=[]
    for a in sorted(adj):
        if a in done:continue
        todo=[a];seen=set()
        while todo:
            b=todo.pop()
            if b in seen:continue
            seen.add(b);todo.extend(adj[b]-seen)
        done|=seen
        sizes.append({'single_mutants':len(seen),'groups':sum(a in seen and b in seen for a,b in pairs)})
    return sorted(sizes,key=lambda x:x['groups'],reverse=True)

def genotype_audit(work,out):
    workbook=work/'programme_sources/gfp_supplementary_data.xlsx'
    rs=xlsx_rows(workbook,7);labels={};rows={};presence={}
    for rn,d in rs[2:]:
        key=mutations(d['A']['text'])
        if key in labels:raise ValueError('Duplicate canonical S7 genotype '+mlabel(key))
        labels[key]=rn;rows[rn]=key
        presence[key]={c for c in COLS if c in d and d[c]['kind']=='n'}
    refs={}; peptides={};cds_sha={}
    for g in ['GFP','URA3']:
        seq=str(SeqIO.read(work/'primary_biology'/('Ref_'+g+'.fasta'),'fasta').seq).upper()
        end=3*(len(str(Seq(seq[:len(seq)//3*3]).translate(to_stop=True)))+1)
        refs[g]=seq[:end];peptides[g]=str(Seq(refs[g]).translate());cds_sha[g]=hashlib.sha256(refs[g].encode()).hexdigest()
    s5=defaultdict(list)
    for rn,d in xlsx_rows(workbook,5):
        if rn>=10 and 'A' in d and 'B' in d:
            s5[d['A']['text']].append(mutations(d['B']['text']))
    coverage=[];groups=[];conditions={}
    for c,(g,prom,med) in COLS.items():
        candidates={k:rn for k,rn in labels.items() if c in presence[k]}
        valid={};reason=Counter();sequence_hash={}
        for key,rn in candidates.items():
            if any(p>len(refs[g]) for a,p,b in key):reason['outside_CDS']+=1;continue
            if any(refs[g][p-1]!=a for a,p,b in key):reason['reference_allele_mismatch']+=1;continue
            seq=list(refs[g])
            for a,p,b in key:seq[p-1]=b
            seq=''.join(seq)
            if str(Seq(seq).translate())!=peptides[g]:reason['protein_changed']+=1;continue
            valid[key]=rn;sequence_hash[key]=hashlib.sha256(seq.encode()).hexdigest()
        pairs=[];position_pairs=set();n_incomplete=0;selected=[]
        for key,rn in valid.items():
            if len(key)!=2:continue
            a,b=key
            if (a,) not in valid or (b,) not in valid:n_incomplete+=1;continue
            la,lb=mlabel((a,)),mlabel((b,));pairs.append((la,lb));position_pairs.add(tuple(sorted([a[1],b[1]])))
            partition='development' if max(a[1],b[1])<=500 else ('validation' if min(a[1],b[1])>550 else 'quarantine')
            group={'reporter':g,'promoter':prom,'medium':med,'column':c,'group':mlabel(key),'single_A':la,'single_B':lb,
                   'WT_relative_fitness_by_definition':1,'row_A':valid[(a,)],'row_B':valid[(b,)],'row_AB':rn,
                   'position_A':a[1],'position_B':b[1],'same_codon':(a[1]-1)//3==(b[1]-1)//3,
                   'stop_codon_changed':max(a[1],b[1])>len(refs[g])-3,'CDS_nt':len(refs[g]),
                   'protein_identical':True,'CDS_WT_sha256':cds_sha[g], 'CDS_A_sha256':sequence_hash[(a,)],
                   'CDS_B_sha256':sequence_hash[(b,)],'CDS_AB_sha256':sequence_hash[key],'planned_partition_J_only':partition if c=='J' else 'untouched'}
            groups.append(group);selected.append(group)
        source_key=f'{prom}-{g} in {med}';source_list=s5[source_key];source_set=set(source_list)
        comp=component_sizes(pairs)
        row={'reporter':g,'promoter':prom,'medium':med,'column':c,'selected_genotypes':len(candidates),
             'selected_single_mutants':sum(len(k)==1 for k in candidates),'selected_doubles':sum(len(k)==2 for k in candidates),
             'synonymous_singles':sum(len(k)==1 for k in valid),'synonymous_doubles':sum(len(k)==2 for k in valid),
             'complete_synonymous_groups':len(pairs),'participating_unique_singles':len(set(x for ab in pairs for x in ab)),
             'unique_doubles':len(pairs),'unique_position_pairs':len(position_pairs),'synonymous_doubles_missing_single':n_incomplete,
             'outside_CDS_exclusions':reason['outside_CDS'],'allele_mismatch_exclusions':reason['reference_allele_mismatch'],
             'protein_change_exclusions':reason['protein_changed'],'S5_rows':len(source_list),'S5_unique_genotypes':len(source_set),
             'S5_duplicate_rows':len(source_list)-len(source_set),'S5_unique_not_in_S7':len(source_set-set(candidates)),
             'S7_not_in_S5_unique':len(set(candidates)-source_set),'dependency_components':len(comp),
             'largest_component_groups':comp[0]['groups'] if comp else 0,'shared_normalized_WT':True}
        coverage.append(row);conditions[c]={'coverage':row,'component_sizes':comp}
    write_csv(out/'synonymous_group_mapping.csv',groups);write_csv(out/'condition_coverage.csv',coverage)
    dump(out/'genotype_audit.json',{'scope':'No numerical fitness outcomes extracted. Numeric-cell presence only.',
         'source_sha256':sha(workbook),'unique_S7_labels':len(labels),'reference_CDS_sha256':cds_sha,'conditions':conditions,
         'cutoff_identifiability':'S5 is the reported measured-genotype list, not an independently documented pre-cutoff universe. Its duplicates are not replicates. Missingness cannot be attributed individually to the cutoff rather than read coverage/mapping. Lost complete-group counts are not identified.'})
    return coverage

def fnum(x):
    try:
        y=float(x)
        return y if math.isfinite(y) else None
    except (TypeError,ValueError):return None

def structural_features(work,out):
    p=work/'primary_biology';rec=SeqIO.read(p/'assets/NC_012920.1.gb','genbank');ref=str(rec.seq).upper()
    functional_ref=str(SeqIO.read(p/'functional_mito_reference.fasta','fasta').seq).upper()
    if functional_ref!=ref:raise ValueError('Functional reference differs from NC_012920.1')
    if [i for i,b in enumerate(ref) if b not in 'ACGT']!=[3106]:raise ValueError('Unexpected reference ambiguity')
    cds={f.qualifiers['gene'][0]:f for f in rec.features if f.type=='CDS'}
    profiles={};align=[]
    with tarfile.open(p/'assets/mito_processed.tar') as t:
        for m in t.getmembers():
            if not m.isfile():continue
            r=list(csv.DictReader(io.StringIO(gzip.decompress(t.extractfile(m).read()).decode()),delimiter='\t'))
            bases=''.join(x['Base'] for x in r).upper();orient='L' if 'Light' in m.name else 'H'
            target=ref.replace('N','') if orient=='H' else str(Seq(ref.replace('N','')).reverse_complement())
            if bases!=target[:len(bases)]:raise ValueError('Population reference mismatch '+m.name)
            name=re.search(r'(WT[12]|KO[12])_',m.name).group(1)
            profiles[(name,orient)]={int(x['Position'])-1:x for x in r}
            align.append({'source':m.name,'positions':len(bases),'reference':'NC_012920.1 minus N3107; reverse-complement for L', 'sequence_mismatches':0})
    clusters={};cluster_map=[]
    with ZipFile(p/'assets/mitoDMS-repo.zip') as z:
        for cond in ['WT','KO']:
            fn=f'mitoDMS-repo/data/Clustering/{cond}_ClusteringDF_Normalized.csv'
            by=defaultdict(list)
            for r in csv.DictReader(io.StringIO(z.read(fn).decode())):by[r['Gene']].append(r)
            for label,rs in by.items():
                unit='CYTB' if label=='NA' else label
                if unit not in UNITS:raise ValueError('Unmapped cluster gene '+unit)
                s=''.join(x['Base'] for x in rs).upper().replace('U','T');reference=ref if unit!='ND6' else str(Seq(ref).reverse_complement())
                anchors=Counter(reference.find(s[i:i+30])-i for i in range(0,len(s)-30,20) if reference.find(s[i:i+30])>=0)
                if not anchors:raise ValueError('No sequence anchor '+unit)
                origin=anchors.most_common(1)[0][0];mismatch=sum(a!=b for a,b in zip(s,reference[origin:origin+len(s)]))
                if mismatch/len(s)>.01:raise ValueError('Cluster sequence agreement <99% '+unit)
                clusters[(cond,unit)]={origin+i:r for i,r in enumerate(rs)}
                cluster_map.append({'condition':cond,'source_gene_label':label,'unit':unit,'sequence_origin0':origin,
                                    'orientation':'reverse_complement' if unit=='ND6' else 'forward','bases':len(s),'mismatches':mismatch,'anchor_support':anchors[origin]})
    def source_index(g,orient):
        r=ref if orient=='H' else str(Seq(ref).reverse_complement())
        return g-(1 if g>r.index('N') else 0)
    def coords(g):
        f=cds[g];strand=f.location.strand
        start=int(f.location.start) if strand==1 else len(ref)-int(f.location.end)
        return start,len(f.extract(rec.seq)),'H' if strand==1 else 'L'
    def bulk(g,window):
        start,length,orient=coords(g);r=ref if orient=='H' else str(Seq(ref).reverse_complement())
        use=range(start,start+(length if window is None else window));pos=[i for i in use if r[i] in 'AC']
        values={}
        for cond in ['WT','KO']:
            vals=[]
            for rep in [1,2]:
                rows=profiles[(cond+str(rep),orient)]
                xs=[fnum(rows[source_index(i,orient)]['Mismatches + Deletions']) for i in pos]
                # The source heading is explicitly verified below; never substitute a different rate.
                if any(x is None for x in xs):raise ValueError('Missing population DMS value')
                vals.append(float(np.mean(xs)))
            values[cond]=float(np.mean(vals))
        return values['KO']-values['WT']
    # Header spelling varies only in spacing; detect the exact named combined mutation field.
    all_headers=next(iter(next(iter(profiles.values())).values())).keys()
    rate_key=next((k for k in all_headers if 'Mismatches' in k and 'Deletions' in k),None)
    if rate_key is None:raise ValueError('Combined mutation frequency absent')
    if rate_key!='Mismatches + Deletions':
        for prof in profiles.values():
            for row in prof.values():row['Mismatches + Deletions']=row[rate_key]
    rows_out=[];coverage=[]
    for window in [15,30,45]:
        for unit,genes in UNITS.items():
            variances=[];valid_unit=True;details=[]
            for gene in genes:
                start,length,orient=coords(gene);w=clusters[('WT',unit)];k=clusters[('KO',unit)]
                ac=[];usable=[];v={'WT':[],'KO':[]}
                for i in range(start,start+window):
                    rw=w.get(i);rk=k.get(i)
                    if rw is None or rk is None:continue
                    if rw['Base'].upper()!=rk['Base'].upper() or rw['Base'].upper() not in 'AC':continue
                    ac.append(i);temp={}
                    for cond,r in [('WT',rw),('KO',rk)]:
                        nums=[fnum(r[x]) for x in ['Cluster_1','Cluster_2','Prop_Cluster1','Prop_Cluster2']]
                        if any(x is None or x<0 for x in nums):break
                        a,b,pa,pb=nums
                        if pa+pb<=0:break
                        pi=pa/(pa+pb);temp[cond]=pi*(1-pi)*(a-b)**2
                    if len(temp)==2:
                        usable.append(i)
                        for cond in v:v[cond].append(temp[cond])
                ok=len(usable)>=8 and len(usable)>=.5*len(ac)
                valid_unit &= ok
                details.append({'gene':gene,'eligible_shared_AC':len(ac),'usable_AC':len(usable),'passed':ok})
                if ok:variances.append(float(np.mean(v['KO'])-np.mean(v['WT'])))
            coverage.append({'unit':unit,'window':window,'passed':valid_unit,'starts':details})
            if valid_unit:
                rows_out.append({'unit':unit,'window':window,'delta_initiation_DMS':float(np.mean([bulk(gene,window) for gene in genes])),
                                 'delta_CDS_DMS':float(np.mean([bulk(gene,None) for gene in genes])),
                                 'delta_cluster_variance':float(np.mean(variances))})
    write_csv(out/'mitochondrial_structural_features.csv',rows_out)
    dump(out/'structural_gate.json',{'population_sequence_checks':align,'functional_reference_exactly_matches_native_reference':True,
         'cluster_sequence_map':cluster_map,'coverage':coverage,'primary_eligible_units':[r['unit'] for r in rows_out if r['window']==30],
         'functional_outcomes_read':False,'interpretation':'Local normalized probing-mixture variance, not thermodynamic conformer occupancy.'})
    return rows_out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--work',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    cov=genotype_audit(a.work,a.out);features=structural_features(a.work,a.out)
    dump(a.out/'summary.json',{'A_complete_groups_by_condition':cov,'B_primary_units':sum(x['window']==30 for x in features),
         'numerical_functional_outcomes_read':False})
    print(json.dumps({'A_conditions':len(cov),'A_total_condition_groups':sum(r['complete_synonymous_groups'] for r in cov),'B_primary_units':sum(x['window']==30 for x in features),'functional_outcomes_read':False}))
if __name__=='__main__':main()
