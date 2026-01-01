"""Frozen two-phase Gill measurement-reproducibility workflow.

Train interprets A/B/C/E/F/G for selected keys only. Validate verifies the pinned
training receipt and its files before interpreting A/D/G. No authors' code runs.
The caller supplies portable protocol/inventory/source paths. No Git required.
"""
from __future__ import annotations
import argparse, gzip, hashlib, io, json, math, platform, re, sys
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from itertools import product, combinations, combinations_with_replacement
from pathlib import Path
import xml.etree.ElementTree as ET
from zipfile import ZipFile
import numpy as np

NS={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
TAG='{'+NS['x']+'}'
HEADERS={'A':'Sequence','B':'CDS_1 Barcode Count','C':'CDS_2 Barcode Count',
         'D':'CDS_3 Barcode Count','E':'CDS_8 Barcode Count','F':'CDS_9 Barcode Count',
         'G':'CDS_10 Barcode Count'}
TRAIN=('B','C','E','F','G');TEST=('D','G')
SEGMENTS=((0,1),(1,4),(4,7),(7,8))
SCENARIOS=list(product(range(3),repeat=2))
SCENARIO_PAIRS=list(combinations_with_replacement(range(9),2))
LIBRARY_PAIRS=list(combinations(range(3),2))
NUMERICAL_FACTOR=128.0

def require(condition,message):
    if not condition:raise ValueError(message)

def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def raw_json(x):return (json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def write_json(p,x):Path(p).write_bytes(raw_json(x))
def write_gz(p,x):Path(p).write_bytes(gzip.compress(raw_json(x),mtime=0))
def read_json(p):return json.loads(Path(p).read_bytes())
def read_gz(p):return json.loads(gzip.decompress(Path(p).read_bytes()))
def code_hash():return digest(Path(__file__))

def _cell(cell,shared):
    if cell.find('x:f',NS) is not None:raise ValueError('Formula input cell')
    typ=cell.attrib.get('t')
    require(typ not in ('b','e','d'),'Invalid input cell type')
    if typ=='inlineStr':
        node=cell.find('x:is',NS);return None if node is None else ''.join(node.itertext())
    v=cell.find('x:v',NS);text=None if v is None else v.text
    if typ=='s' and text is not None:
        try:return ''.join(shared[int(text)].itertext())
        except (ValueError,IndexError) as e:raise ValueError('Bad selected shared-string index') from e
    return text

def integer(text):
    try:
        x=Decimal(text)
        require(x.is_finite() and x>=0 and x==x.to_integral_value(),'Invalid count')
        return int(x)
    except (InvalidOperation,TypeError,ValueError) as e:raise ValueError('Count must be finite nonnegative integer') from e

def read_selected(path,keys,columns,expected_rows=None):
    """Inspect A first; do not decode any unselected row or forbidden value."""
    wanted=set(keys);columns=tuple(columns);allowed={'A',*columns};result={};decoded=defaultdict(int)
    require(len(wanted)==len(keys),'Duplicate frozen keys')
    with ZipFile(path) as z:
        shared=list(ET.fromstring(z.read('xl/sharedStrings.xml'))) if 'xl/sharedStrings.xml' in z.namelist() else []
        sheets=[n for n in z.namelist() if re.fullmatch(r'xl/worksheets/sheet\d+\.xml',n)]
        require(len(sheets)==1,'Expected exactly one source worksheet')
        with z.open(sheets[0]) as stream:
            for _,row in ET.iterparse(stream,events=('end',)):
                if row.tag!=TAG+'row':continue
                cells={}
                for c in row:
                    col=re.sub(r'\d','',c.attrib.get('r',''))
                    require(col not in cells,'Duplicate cell column')
                    cells[col]=c
                rn=int(row.attrib['r'])
                if rn==1:
                    for col in allowed:
                        require(col in cells,'Missing selected header')
                        require(_cell(cells[col],shared)==HEADERS[col],'Selected header mismatch '+col)
                        decoded[col]+=1
                else:
                    require('A' in cells,'Missing sequence identity')
                    key=_cell(cells['A'],shared);decoded['A']+=1
                    if key in wanted:
                        require(key not in result,'Duplicate selected genotype')
                        if expected_rows is not None:require(rn==expected_rows[key],'Frozen row identity mismatch')
                        values={}
                        for col in columns:
                            require(col in cells,'Missing selected count '+col)
                            values[col]=integer(_cell(cells[col],shared));decoded[col]+=1
                        result[key]=values
                row.clear()
    require(set(result)==wanted,'Missing selected genotype rows')
    return result,dict(decoded)

def design(item):
    keys=item['observed_keys'];n=len(keys)
    require(keys==sorted(set(keys)) and n==item['expected_vertices'],'Frozen membership invalid')
    require(item['complete'] and len(item['expected_component_sizes'])==1,'Incomplete/disconnected frozen class')
    levels=[sorted({k[a:b] for k in keys}) for a,b in SEGMENTS]
    require([''.join(v) for v in product(*levels)]==keys,'Incomplete codon Cartesian product')
    p=1+sum(len(x)-1 for x in levels)
    require(p==item['codon_main_effect_rank'] and n-p==item['codon_residual_dimensions'] and n-p>=3,'Frozen rank mismatch')
    groups=[]
    for (a,b),dom in zip(SEGMENTS,levels):
        groups.append([np.array([i for i,k in enumerate(keys) if k[a:b]==v],dtype=int) for v in dom])
    adjacency=[[] for _ in keys]
    for i,x in enumerate(keys):
        for j in range(i):
            if sum(a!=b for a,b in zip(x,keys[j]))==1:adjacency[i].append(j);adjacency[j].append(i)
    require(sum(map(len,adjacency))//2==item['expected_edges'],'Frozen edge count mismatch')
    return {'keys':keys,'groups':groups,'adjacency':adjacency,'n':n,'p':p,'q':n-p}

def project(f,d):
    f=np.asarray(f,dtype=float);require(f.shape==(d['n'],) and np.all(np.isfinite(f)),'Invalid field')
    mean=math.fsum(map(float,f))/len(f);p=np.full(len(f),-3*mean)
    for groups in d['groups']:
        for idx in groups:p[idx]+=math.fsum(float(f[i]) for i in idx)/len(idx)
    q=f-p
    # Fixed numerical zero rule, never a measurement or statistical cutoff.
    tol=NUMERICAL_FACTOR*np.finfo(float).eps*len(f)*max(1.0,float(np.max(np.abs(f))))
    q[np.abs(q)<=tol]=0.0
    return p,q

def centered(x):
    x=np.asarray(x,dtype=float);return x-math.fsum(map(float,x))/len(x)
def dot(a,b):return math.fsum(float(x)*float(y) for x,y in zip(a,b))
def cosine(a,b):
    aa=dot(a,a);bb=dot(b,b)
    if aa==0 or bb==0:return None
    r=dot(a,b)/math.sqrt(aa*bb)
    require(abs(r)<=1+1e-12,'Cosine outside numerical bound')
    return max(-1.0,min(1.0,r))
def ranks(x):
    x=list(x);order=sorted(range(len(x)),key=lambda i:(x[i],i));out=np.empty(len(x));i=0
    while i<len(order):
        j=i+1
        while j<len(order) and x[order[j]]==x[order[i]]:j+=1
        for k in order[i:j]:out[k]=(i+j-1)/2
        i=j
    return out
def spearman(a,b):return cosine(centered(ranks(a)),centered(ranks(b)))

def maximum_barrier(f,adj):
    """Exact all-pair maximum via superlevel component merges, ties allowed."""
    f=list(map(float,f));n=len(f);parent=list(range(n));peak=list(range(n));active=[False]*n;best=0.0;pair=None
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for v in sorted(range(n),key=lambda i:(-f[i],i)):
        active[v]=True
        for u in adj[v]:
            if not active[u]:continue
            a,b=find(v),find(u)
            if a==b:continue
            x,y=peak[a],peak[b];bar=min(f[x],f[y])-f[v]
            candidate=tuple(sorted((x,y)))
            if bar>best or (bar==best and bar>0 and (pair is None or candidate<pair)):best=bar;pair=candidate
            parent[b]=a;peak[a]=min((x,y),key=lambda i:(-f[i],i))
    require(len({find(i) for i in range(n)})==1,'Barrier graph disconnected')
    precision_tolerance=float(NUMERICAL_FACTOR*np.finfo(float).eps*n*max(1.0,max(map(abs,f))))
    return {'maximum':max(0.0,best),'precision_tolerance':precision_tolerance,'numerically_resolved_positive':bool(best>precision_tolerance),'unresolved_tiny_positive':bool(0<best<=precision_tolerance),'arg_pair':None if pair is None else list(pair),'endpoint_pairs_covered':n*(n-1)//2}

def load_inputs(protocol_path,inventory_path,expected_protocol_hash):
    require(digest(protocol_path)==expected_protocol_hash,'Protocol hash mismatch')
    protocol=read_json(protocol_path);require(protocol['id']=='GILL_MULTICLASS_REPRODUCIBILITY_V1' and protocol['status']=='frozen_before_nonreference_RNA_outcomes','Unrecognized frozen protocol');raw=gzip.decompress(Path(inventory_path).read_bytes())
    require(hashlib.sha256(raw).hexdigest()==protocol['inventory_sha256'],'Inventory content hash mismatch')
    inv=json.loads(raw);lookup={x['class_id']:x for x in inv['classes']}
    require(len(lookup)==len(inv['classes']),'Duplicate inventory class')
    classes=[];parts=protocol['class_ids_by_partition'];ids=parts['development']+parts['validation']
    require(len(set(ids))==len(ids) and not(set(ids)&set(parts['reference_stratum_quarantine'])),'Partition overlap')
    for part in ['development','validation']:
        for cid in parts[part]:
            require(cid in lookup,'Missing frozen class');item=lookup[cid]
            require(item['candidate_reproducibility_eligible'] and item['candidate_partition']==part and not item['reference_protein'],'Ineligible selected class')
            classes.append((part,item,design(item)))
        require(len(parts[part])==protocol[('primary' if part=='validation' else 'development')+'_class_count'],'Frozen class denominator mismatch')
        require(sum(d['n'] for p,_,d in classes if p==part)==protocol[('primary' if part=='validation' else 'development')+'_vertices'],'Frozen vertex denominator mismatch')
    keys=[k for _,_,d in classes for k in d['keys']];require(len(set(keys))==len(keys),'Overlapping genotype classes')
    require(protocol['pseudocounts']==[.25,.5,1.0] and protocol['primary_pseudocount']==.5,'Unexpected pseudocount settings')
    require(protocol['training_RNA_columns']==['B','C'] and protocol['training_DNA_columns']==['E','F'] and protocol['validation_RNA_column']=='D' and protocol['validation_DNA_column']=='G','Unexpected library settings')
    hashes={'protocol_sha256':digest(protocol_path),'inventory_sha256':protocol['inventory_sha256'],'inventory_gzip_sha256':digest(inventory_path),'code_sha256':code_hash()}
    return protocol,classes,hashes

def logs(counts,keys,cols,c):
    x=np.array([[math.log2(counts[k][col]+c) for k in keys] for col in cols],dtype=float)
    require(np.all(np.isfinite(x)),'Nonfinite transformed count')
    return x

def training_fields(counts,classes):
    out={}
    for part,item,d in classes:
        fields={}
        for c in [.25,.5,1.0]:
            lr=logs(counts,d['keys'],('B','C'),c);ld=logs(counts,d['keys'],('E','F','G'),c)
            r=np.array([project(v,d)[1] for v in lr]);z=np.array([project(v,d)[1] for v in ld])
            t=(lr[0]+lr[1]-ld[0]-ld[1])/2;pt,qt=project(t,d);u=(r[0]+r[1]-z[0]-z[1])/2
            require(np.max(np.abs(qt-u))<=1e-10*max(1.,float(np.max(np.abs(t)))),'Projection linearity discrepancy')
            fields[str(c)]={'T':t.tolist(),'PT':pt.tolist(),'QT':u.tolist(),'RNA_raw_counts':[[counts[k][col] for k in d['keys']] for col in ('B','C')],'DNA_raw_counts':[[counts[k][col] for k in d['keys']] for col in ('E','F','G')],'RNA_logs':lr.tolist(),'DNA_logs':ld.tolist(),'RNA_residuals':r.tolist(),'DNA_residuals':z.tolist(),'training_barrier':maximum_barrier(t,d['adjacency']),'training_codon_barrier':maximum_barrier(pt,d['adjacency'])}
        out[item['class_id']]={'partition':part,'keys':d['keys'],'fields':fields}
    return out

def fresh(path):
    path=Path(path);require(not path.exists(),'Output directory already exists');path.mkdir(parents=True);return path
def record_inputs(workbook,protocol,hashes):
    require(digest(workbook)==protocol['source_workbook_sha256'],'Workbook hash mismatch before protected read')
    return dict(hashes,workbook_sha256=protocol['source_workbook_sha256'])

def train(args):
    protocol,classes,hashes=load_inputs(args.protocol,args.inventory,args.expected_protocol_sha256)
    inp=record_inputs(args.workbook,protocol,hashes);keys=sorted(k for _,_,d in classes for k in d['keys'])
    expected={k:r for _,g,_ in classes for k,r in zip(g['observed_keys'],g['observed_workbook_rows'])}
    require(not Path(args.output).exists(),'Output directory already exists')
    counts,exposure=read_selected(args.workbook,keys,TRAIN,expected);pred=training_fields(counts,classes)
    require(digest(args.workbook)==inp['workbook_sha256'],'Workbook changed during train')
    out=fresh(args.output);write_gz(out/'counts.json.gz',counts);write_gz(out/'predictions.json.gz',pred)
    summary={'phase':'train','classes':len(classes),'vertices':len(keys),'decoded_columns':exposure,'D_value_or_header_decoded':False,'settings':{'zero_rule_factor':NUMERICAL_FACTOR,'pseudocounts':[.25,.5,1.0]},'scope':'conditional processed-library reproduction; no biological calibration'}
    write_json(out/'training_summary.json',summary)
    receipt={'phase':'train','inputs':inp,'outputs':{f:digest(out/f) for f in ['counts.json.gz','predictions.json.gz','training_summary.json']},'runtime':{'python':platform.python_version(),'numpy':np.__version__},'D_opened':False}
    write_json(out/'receipt.json',receipt);return receipt

def check_training(args,protocol,classes,hashes):
    base=Path(args.training);rp=base/'receipt.json'
    require(digest(rp)==args.expected_training_receipt_sha256,'Pinned training receipt mismatch')
    receipt=read_json(rp);expected=record_inputs(args.workbook,protocol,hashes)
    require(receipt['phase']=='train' and receipt['D_opened'] is False and receipt['inputs']==expected,'Training identity mismatch')
    names={'counts.json.gz','predictions.json.gz','training_summary.json'}
    require(set(receipt['outputs'])==names,'Unexpected training receipt output set')
    for f in names:require(digest(base/f)==receipt['outputs'][f],'Training artifact hash mismatch '+f)
    counts=read_gz(base/'counts.json.gz');pred=read_gz(base/'predictions.json.gz')
    require(set(counts)=={k for _,_,d in classes for k in d['keys']},'Stored count membership changed')
    for row in counts.values():
        require(set(row)==set(TRAIN),'Stored count columns')
        for v in row.values():require(isinstance(v,int) and not isinstance(v,bool) and v>=0,'Stored invalid count')
    recomputed=training_fields(counts,classes)
    require(raw_json(recomputed)==raw_json(pred),'Frozen predictions disagree with saved training counts/runtime')
    return receipt,counts,pred

def quad_terms(a,b,c,d):
    terms=[dot(a,c),dot(b,d),-dot(a,d),-dot(b,c)];total=dot(a-b,c-d)
    err=abs(math.fsum(terms)-total);scale=max(1.0,sum(map(abs,terms)),abs(total))
    require(err<=1e-11*scale,'Quadratic decomposition numerical failure')
    return terms+[total],err

def class_validation(item,d,f,ld3,raw_rna3):
    lr=np.vstack([np.array(f['RNA_logs']),ld3]);dna=np.array(f['DNA_logs'])
    rawrna=f['RNA_raw_counts']+[raw_rna3];rawdna=f['DNA_raw_counts']
    r=np.vstack([np.array(f['RNA_residuals']),project(ld3,d)[1]]);z=np.array(f['DNA_residuals'])
    t=np.array(f['T']);u=np.array(f['QT']);vraw=lr[2]-dna[2];pv,qv=project(vraw,d);v=r[2]-z[2]
    require(np.max(np.abs(qv-v))<=1e-10*max(1.,float(np.max(np.abs(vraw)))),'Held-out projection linearity discrepancy')
    terms=[];errors=[]
    for i,j in SCENARIO_PAIRS:
        a,b=SCENARIOS[i];c,e=SCENARIOS[j];q,err=quad_terms(r[a],z[b],r[c],z[e]);terms.append(q);errors.append(err)
    primary,pe=quad_terms((r[0]+r[1])/2,(z[0]+z[1])/2,r[2],z[2])
    n=d['n'];stats={'cross':dot(u,v)/n,'train_energy':dot(u,u)/n,'validation_energy':dot(v,v)/n,'residual_prediction_mse':dot(u-v,u-v)/n,'zero_residual_mse':dot(v,v)/n,'cosine':cosine(u,v),'whole_centered_cosine':cosine(centered(t),centered(vraw)),'whole_prediction_mse':dot(centered(t)-centered(vraw),centered(t)-centered(vraw))/n,'constant_prediction_mse':dot(centered(vraw),centered(vraw))/n,'codon_prediction_mse':dot(centered(np.array(f['PT']))-centered(vraw),centered(np.array(f['PT']))-centered(vraw))/n}
    barriers=[];additive=[]
    for a,b in SCENARIOS:
        field=lr[a]-dna[b];p,_=project(field,d);barriers.append(maximum_barrier(field,d['adjacency']));additive.append(maximum_barrier(p,d['adjacency']))
    return {'class_id':item['class_id'],'stratum':item['central_amino_acid_stratum_157_158'],'vertices':n,'codon_rank':d['p'],'residual_dimensions':d['q'],'clique_product':item['codon_factors_all_cliques'],'stats':stats,'RNA_residual_gram':[[dot(a,b)/n for b in r] for a in r],'DNA_residual_gram':[[dot(a,b)/n for b in z] for a in z],'RNA_residual_cosines':[cosine(r[a],r[b]) for a,b in LIBRARY_PAIRS],'DNA_residual_cosines':[cosine(z[a],z[b]) for a,b in LIBRARY_PAIRS],'RNA_rank_correlations':[spearman(rawrna[a],rawrna[b]) for a,b in LIBRARY_PAIRS],'DNA_rank_correlations':[spearman(rawdna[a],rawdna[b]) for a,b in LIBRARY_PAIRS],'decomposition_terms':terms,'primary_decomposition_terms':primary,'decomposition_max_abs_error':max(errors+[pe]),'scenario_barriers':barriers,'scenario_codon_barriers':additive,'heldout_barrier':maximum_barrier(vraw,d['adjacency']),'heldout_codon_barrier':maximum_barrier(pv,d['adjacency'])}

def pooled(stats):
    cross=math.fsum(x['cross'] for x in stats);a=math.fsum(x['train_energy'] for x in stats);b=math.fsum(x['validation_energy'] for x in stats)
    rr=None if a==0 or b==0 else cross/math.sqrt(a*b)
    if rr is not None:require(abs(rr)<=1+1e-12,'Pooled cosine outside numerical bound');rr=max(-1.0,min(1.0,rr))
    return {'cosine':rr,'cross':cross,'train_energy':a,'validation_energy':b,'class_count':len(stats)}

def summarize(records,target):
    p=pooled([x['stats'] for x in records]);strata=defaultdict(list)
    for x in records:strata[x['stratum']].append(x)
    influence=[]
    for s,rows in sorted(strata.items()):
        removed=pooled([x['stats'] for x in rows]);rem=pooled([x['stats'] for x in records if x['stratum']!=s])
        influence.append({'stratum':s,'remaining_cosine':rem['cosine'],'train_energy_fraction':None if p['train_energy']==0 else removed['train_energy']/p['train_energy'],'validation_energy_fraction':None if p['validation_energy']==0 else removed['validation_energy']/p['validation_energy']})
    finite=[x['remaining_cosine'] for x in influence if x['remaining_cosine'] is not None]
    fractions=lambda k:[x[k] for x in influence if x[k] is not None]
    aggregate_terms=np.zeros((45,5));primary_terms=np.zeros(5)
    for x in records:aggregate_terms+=np.array(x['decomposition_terms'])/x['vertices'];primary_terms+=np.array(x['primary_decomposition_terms'])/x['vertices']
    losses={k:math.fsum(x['stats'][k] for x in records)/len(records) for k in ['residual_prediction_mse','zero_residual_mse','whole_prediction_mse','constant_prediction_mse','codon_prediction_mse']}
    rgram=np.sum([x['RNA_residual_gram'] for x in records],axis=0);dgram=np.sum([x['DNA_residual_gram'] for x in records],axis=0)
    gramcos=lambda g,a,b:None if g[a,a]==0 or g[b,b]==0 else max(-1.,min(1.,float(g[a,b]/math.sqrt(g[a,a]*g[b,b]))))
    diagonal={a:aggregate_terms[i,4] for i,(a,b) in enumerate(SCENARIO_PAIRS) if a==b}
    scenario_cos=[None if diagonal[a]==0 or diagonal[b]==0 else max(-1.,min(1.,float(aggregate_terms[i,4]/math.sqrt(diagonal[a]*diagonal[b])))) for i,(a,b) in enumerate(SCENARIO_PAIRS)]
    barrier_summary={}
    for architecture,subset in [('clique',[x for x in records if x['clique_product']]),('nonclique',[x for x in records if not x['clique_product']])]:
        barrier_summary[architecture]={'classes':len(subset)}
        for name in ['scenario_barriers','scenario_codon_barriers']:
            barrier_summary[architecture][name]=[{'mean':math.fsum(x[name][i]['maximum'] for x in subset)/len(subset),'maximum':max(x[name][i]['maximum'] for x in subset),'raw_positive_classes':sum(x[name][i]['maximum']>0 for x in subset),'numerically_resolved_positive_classes':sum(x[name][i]['numerically_resolved_positive'] for x in subset),'unresolved_tiny_positive_classes':sum(x[name][i]['unresolved_tiny_positive'] for x in subset)} for i in range(9)] if subset else []
    rank_summary={}
    for name in ['RNA_rank_correlations','DNA_rank_correlations']:
        rank_summary[name]=[]
        for i in range(3):
            vals=[x[name][i] for x in records if x[name][i] is not None]
            rank_summary[name].append({'defined':len(vals),'undefined':len(records)-len(vals),'mean_defined':None if not vals else math.fsum(vals)/len(vals),'range_defined':None if not vals else [min(vals),max(vals)]})
    p.update({'RNA_residual_pooled_cosines':[gramcos(rgram,a,b) for a,b in LIBRARY_PAIRS],'DNA_residual_pooled_cosines':[gramcos(dgram,a,b) for a,b in LIBRARY_PAIRS],'scenario_residual_pooled_cosines':scenario_cos,'rank_agreement_summary':rank_summary,'barriers_by_architecture':barrier_summary,'status':'information_limited' if p['cosine'] is None else 'target_met' if p['cosine']>=target else 'target_failed','target':target,'losses':losses,'leave_one_stratum_out':influence,'leave_one_stratum_out_range':None if not finite else [min(finite),max(finite)],'largest_stratum_training_energy_fraction':max(fractions('train_energy_fraction'),default=None),'largest_stratum_validation_energy_fraction':max(fractions('validation_energy_fraction'),default=None),'decomposition_aggregate_terms':aggregate_terms.tolist(),'primary_decomposition_aggregate_terms':primary_terms.tolist(),'zero_energy_classes':{'training':sum(x['stats']['train_energy']==0 for x in records),'validation':sum(x['stats']['validation_energy']==0 for x in records)},'clique_classes':sum(x['clique_product'] for x in records),'nonclique_classes':sum(not x['clique_product'] for x in records)})
    return p

def validate(args):
    protocol,classes,hashes=load_inputs(args.protocol,args.inventory,args.expected_protocol_sha256)
    require(not Path(args.output).exists(),'Output directory already exists')
    receipt,counts,pred=check_training(args,protocol,classes,hashes) # All before D.
    keys=sorted(counts);expected={k:r for _,g,_ in classes for k,r in zip(g['observed_keys'],g['observed_workbook_rows'])}
    test,exposure=read_selected(args.workbook,keys,TEST,expected)
    for k in keys:require(test[k]['G']==counts[k]['G'],'Historical DNA column changed')
    records={};summaries={}
    for c in [.25,.5,1.0]:
        split={'development':[],'validation':[]}
        for part,item,d in classes:
            lr3=logs(test,d['keys'],('D',),c)[0];q=class_validation(item,d,pred[item['class_id']]['fields'][str(c)],lr3,[test[k]['D'] for k in d['keys']]);split[part].append(q)
        records[str(c)]=split;summaries[str(c)]={part:summarize(rs,protocol['primary_descriptive_reproducibility_target']) for part,rs in split.items()}
    require(digest(args.workbook)==protocol['source_workbook_sha256'],'Workbook changed during validation')
    out=fresh(args.output);write_gz(out/'heldout_counts.json.gz',test);write_gz(out/'class_metrics.json.gz',records)
    result={'phase':'validate','primary_pseudocount':.5,'primary_partition':'validation','primary_result':summaries['0.5']['validation'],'all_summaries':summaries,'scenario_order':SCENARIOS,'decomposition_pair_order':SCENARIO_PAIRS,'library_pair_order':LIBRARY_PAIRS,'decomposition_term_order':['RNA_RNA','DNA_DNA','minus_RNA_DNA','minus_DNA_RNA','total'],'exposure_decoded_columns':exposure,'uncertainty':protocol['uncertainty'],'response_scale_limit':protocol['response_scale_limit'],'centered_loss_interpretation':'Whole-field constant/codon/T losses compare shape after separately centering T and V; validation centering uses held-out values to define this diagnostic, so these are not out-of-sample level predictions.','barrier_precision_rule':'Raw maximum is retained. Numerical resolution uses128*machine_epsilon*n*max(1,maxabs(field)); this is a conservative machine-precision diagnostic, not measurement uncertainty or a statistical threshold.','numerical_zero_rule':{'factor':NUMERICAL_FACTOR,'formula':'factor*machine_epsilon*n*max(1,maxabs(field))'},'retained_failed_reference_pilot':True}
    write_json(out/'validation.json',result)
    final={'phase':'validate','inputs':dict(hashes,workbook_sha256=protocol['source_workbook_sha256'],training_receipt_sha256=args.expected_training_receipt_sha256),'outputs':{f:digest(out/f) for f in ['heldout_counts.json.gz','class_metrics.json.gz','validation.json']},'D_opened_after_training_gate':True,'runtime':{'python':platform.python_version(),'numpy':np.__version__}}
    write_json(out/'receipt.json',final);return final

def main():
    if not __debug__:raise RuntimeError('Optimized Python disables required verification; use ordinary Python')
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['train','validate'])
    for name in ['protocol','inventory','workbook','output']:p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--expected-protocol-sha256',required=True);p.add_argument('--training',type=Path);p.add_argument('--expected-training-receipt-sha256')
    a=p.parse_args()
    if a.phase=='validate':require(a.training is not None and a.expected_training_receipt_sha256 is not None,'Validation requires explicitly pinned training receipt')
    result=train(a) if a.phase=='train' else validate(a)
    print(json.dumps({'phase':result['phase'],'output':str(a.output),'receipt_sha256':digest(a.output/'receipt.json')}))

if __name__=='__main__':main()
