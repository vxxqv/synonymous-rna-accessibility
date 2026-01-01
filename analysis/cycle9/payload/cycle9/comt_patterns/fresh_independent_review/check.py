"""Fresh exact COMT pair reconstruction, without analysis/verifier imports."""
import sys
if not __debug__:raise RuntimeError('Optimized Python refused before access')
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from collections import Counter,defaultdict
import hashlib,json,math,statistics
R=Path.cwd();O=Path(__file__).resolve().parent;B=O.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def need(x,m):
    if not x:raise ValueError(m)
def number(x):need(type(x) in (int,float) and math.isfinite(x),'invalid stored number');return F(x)
def label(points):
    # A single direction must be strict at every shared-D joint point.
    signs={(int(r>0)-int(r<0),int(t>0)-int(t<0)) for r,t in points}
    return {(1,-1):'reversal_R_up_T_down',(-1,1):'reversal_R_down_T_up',(1,1):'same_direction_a_higher',(-1,-1):'same_direction_b_higher'}.get(next(iter(signs)),'unresolved_with_ties' if any(0 in z for z in signs) else 'unresolved_no_ties') if len(signs)==1 else ('unresolved_with_ties' if any(0 in z for z in signs) else 'unresolved_no_ties')
def cases(z):return [(z['R'][i]-z['D'][j],z['D'][j]-z['F'][k]) for i,j,k in product(range(4),repeat=3)]
def synthetic():
    records=[]
    def ok(z,m):need(z,m);records.append(m)
    z={'R':[F(0)]*4,'D':[F(-1)]*4,'F':[F(0)]*4};p=cases(z)
    ok(len(p)==64 and label(p*2)=='reversal_R_up_T_down','shared-D counterexample and128cases')
    ok(label([(-r,-t) for r,t in p])=='reversal_R_down_T_up','lexical orientation swap preserves reversal')
    ok(label([(F(1),F(1))]*128)=='same_direction_a_higher','same direction not reversal')
    ok(label(p+[(F(-1),F(1))]*64)=='unresolved_no_ties','opposite map orientations do not certify')
    ok(label(p[:-1]+[(F(0),F(-1))])=='unresolved_with_ties','exact zero is not strict')
    ok(label(p[:-1]+[(F(-1,10**30),F(-1))])=='unresolved_no_ties','one arbitrarily small contrary case defeats universal predicate')
    # Primitive per-allele perturbations attain4delta; sharedD is preserved.
    for delta in (F(1,5),F(1,4),F(3,10)):
        # zR falls by2delta and zD rises by2delta at the same group.
        r=F(1)-4*delta;t=F(-1)+2*delta
        ok((r>0)==(delta<F(1,4)),'sharp m/4 boundary '+str(delta))
    aa=[F(1),F(2),F(3),F(4)];bb=[F(0),F(0),F(0),F(0)]
    ok([a-b for a,b in zip(aa,bb)]==[(a+11)-(b+11) for a,b in zip(aa,bb)],'common assay offset cancellation')
    # Source negative-control subtraction can prevent WTdenominator cancellation.
    ratio=lambda Z:(F(4)-Z*F(1,10))/(F(2)-Z*F(1,5))
    ok(ratio(F(1))!=ratio(F(2)),'variant-specific control prevents complete WT cancellation')
    return records
def main():
    out=O/'run_v1';need(not out.exists(),'fresh output required');tests=synthetic()
    freeze=json.loads((B/'execution_freeze.json').read_text());pins={str(B/'execution_freeze.json'):sha(B/'execution_freeze.json')}
    for name,pin in freeze['owned_files'].items():need(sha(B/name)==pin,'changed frozen source '+name);pins[str(B/name)]=pin
    for role,path in freeze['input_paths'].items():need(sha(R/path)==freeze['input_hashes'][role],'changed input '+role);pins[str(R/path)]=sha(R/path)
    need(sha(B/'observed_v1/result.json')=='01136e43574d5d0e0e3b8bb649dc48cd5ac9dbe034ce3b9c6d79096c278b4727','result hash changed');pins[str(B/'observed_v1/result.json')]=sha(B/'observed_v1/result.json')
    contract=json.loads((R/freeze['input_paths']['contract']).read_text());data=json.loads((R/freeze['input_paths']['input']).read_text());inv=json.loads((B/'metadata_inventory.json').read_text())
    targets={x['hgvs_nt']:x for x in contract['targets']};need(len(targets)==125,'target denominator');members={x['sample']:x for x in contract['members'].values()}
    index={};checked_groups=0;median_float_differences=[]
    for mp in ('code','geo'):
        index[mp]={r['hgvs_nt']:r for r in data[mp]};need(len(data[mp])==len(index[mp])==125 and set(index[mp])==set(targets),'source identity coverage')
        for ident,row in index[mp].items():
            for field in ('ROI','site','native_codon_site','is_SNV','edit_key'):need(row[field]==targets[ident][field],'source metadata disagreement')
            for assay in ('R','D','F'):
                gs=row['assays'][assay]['groups'];need(len(gs)==4 and {g['bio'] for g in gs}=={1,2,3,4},'group completeness')
                for g in gs:
                    samples,values=g['samples'],g['values'];need(len(samples)==len(values)==len(set(samples)),'sample value alignment')
                    mask=[x is not None for x in values];need(sum(mask)==g['technical_finite'] and sum(mask)>0,'finite support mismatch')
                    finite=[x for x in values if x is not None];[number(x) for x in finite];number(g['median'])
                    median_float_differences.append(abs(statistics.median(finite)-g['median']))
                    for s in samples:need((members[s]['trait'],members[s]['ROI'],members[s][mp+'_bio'])==(assay,row['ROI'],g['bio']),'contract group mapping mismatch')
                    checked_groups+=1
    expected=[];pair_checks=[];allpoints={}
    for a,b in combinations(sorted(targets,key=lambda s:s.encode('utf-8')),2):
        if targets[a]['native_codon_site']!=targets[b]['native_codon_site']:continue
        meta={'a':a,'b':b,'site':targets[a]['native_codon_site'],'ROI':targets[a]['ROI'],'both_SNV':targets[a]['is_SNV'] and targets[b]['is_SNV']}
        permap={};points=[]
        for mp in ('code','geo'):
            z={};sources={}
            for assay in ('R','D','F'):
                aa={g['bio']:g for g in index[mp][a]['assays'][assay]['groups']};bb={g['bio']:g for g in index[mp][b]['assays'][assay]['groups']};z[assay]=[];sources[assay]=[]
                for group in range(1,5):
                    x,y=aa[group],bb[group];need(x['samples']==y['samples'],'pair sample list mismatch');need([v is not None for v in x['values']]==[v is not None for v in y['values']],'pair contribution-mask mismatch')
                    z[assay].append(number(x['median'])-number(y['median']));sources[assay].append((x['samples'],[v is not None for v in x['values']]))
            pp=cases(z);permap[mp]={'z':z,'sources':sources,'points':pp,'category':label(pp)};points+=pp
        category=label(points);bounds={'R_lower':min(x[0] for x in points),'R_upper':max(x[0] for x in points),'T_lower':min(x[1] for x in points),'T_upper':max(x[1] for x in points)}
        margin=min(abs(x) for p in points for x in p) if category.startswith('reversal') else F(0)
        expected.append({**meta,'category':category,'complete':True,'combined_bounds':{k:str(x) for k,x in bounds.items()},'unstructured_log_radius_exclusive':str(margin/4),'map_categories_differ':len({x['category'] for x in permap.values()})>1})
        allpoints[a,b]=permap
    need(len(expected)==134,'pair count')
    need({(r['a'],r['b']) for r in expected}=={(r['a'],r['b']) for r in inv['pairs']},'inventory mismatch')
    # Decode answers only after independent exact records are complete.
    answers=json.loads((B/'observed_v1/result.json').read_text());answer={(r['a'],r['b']):r for r in answers['pairs']};need(len(answer)==134,'answer coverage');joint_checks=0;source_masks=0;extremum_checks=0
    for own in expected:
        a,b=own['a'],own['b'];saved=answer[a,b]
        for k,x in own.items():need(saved[k]==x,'pair scalar mismatch '+k+' '+a+' '+b)
        for mp,d in allpoints[a,b].items():
            s=saved['maps'][mp];need(s['complete'] and not s['issues'] and s['category']==d['category'],'permap status mismatch')
            for assay in ('R','D','F'):
                need(s['same_assay_differences'][assay]==[str(x) for x in d['z'][assay]],'same-assay value mismatch')
                for row,(samples,mask) in zip(s['source_groups'][assay],d['sources'][assay]):need(row['samples_a']==row['samples_b']==samples and row['finite_mask_a']==row['finite_mask_b']==mask and row['issue'] is None,'pair source-mask lineage mismatch');source_masks+=1
            need(len(s['joint_cases'])==64,'joint record count')
            for saved_point,(i,j,k),(rr,tt) in zip(s['joint_cases'],product(range(1,5),repeat=3),d['points']):
                need(saved_point=={'R_group':i,'D_group':j,'F_group':k,'DeltaR':str(rr),'DeltaT':str(tt)},'sharedD joint record mismatch');joint_checks+=1
            bounds={'R_lower':min(x[0] for x in d['points']),'R_upper':max(x[0] for x in d['points']),'T_lower':min(x[1] for x in d['points']),'T_upper':max(x[1] for x in d['points'])}
            need(s['bounds']=={k:str(x) for k,x in bounds.items()},'map bounds mismatch')
            for field,bound in bounds.items():
                axis=0 if field.startswith('R') else 1;idx=[i for i,p in enumerate(d['points']) if p[axis]==bound];need(s['extremum_case_indices_zero_based'][field]==idx,'extremum tie provenance mismatch');extremum_checks+=1
            need(s['any_exact_zero']==any(x==0 for p in d['points'] for x in p),'exact zero flag mismatch')
    summaries={}
    for subset in ('all','SNV'):
        rr=expected if subset=='all' else [r for r in expected if r['both_SNV']];vv=list(targets.values()) if subset=='all' else [r for r in targets.values() if r['is_SNV']];by=defaultdict(list)
        for t in vv:by[t['native_codon_site']].append(t)
        site_rows=[]
        for site,ts in sorted(by.items()):
            p=[r for r in rr if r['site']==site];n=sum(r['category'].startswith('reversal') for r in p)
            site_rows.append({'site':site,'variant_count':len(ts),'pair_count':len(p),'complete_pairs':len(p),'uniform_reversals':n,'category_counts':dict(sorted(Counter(r['category'] for r in p).items())),'reversal_fraction':str(F(n,len(p))) if p else None})
        informative=[x for x in site_rows if x['pair_count']];summary={'variants':len(vv),'sites':len(by),'pairs':len(rr),'sites_with_pairs':len(informative),'singleton_sites':len(by)-len(informative),'complete_pairs':len(rr),'joint_case_records':128*len(rr),'category_counts':dict(sorted(Counter(r['category'] for r in rr).items())),'uniform_reversal_pairs':sum(r['category'].startswith('reversal') for r in rr),'sites_with_uniform_reversal':sum(x['uniform_reversals']>0 for x in site_rows),'equal_informative_site_reversal_fraction':str(sum((F(x['reversal_fraction']) for x in informative),F(0))/len(informative)),'map_category_disagreements':sum(r['map_categories_differ'] for r in rr),'pair_counts_by_ROI':dict(sorted(Counter(r['ROI'] for r in rr).items())),'sites_detail':site_rows}
        # JSONobject keys become strings in the saved record.
        need(json.loads(json.dumps(summary))==answers['summaries'][subset],'summary mismatch '+subset);summaries[subset]=summary
    for p,pin in pins.items():need(sha(Path(p))==pin,'input changed during review')
    named=[x for x in expected if x['category'].startswith('reversal')];need(len(named)==1,'named result changed')
    out.mkdir();result={'status':'PASS_FRESH_EXACT_RECONSTRUCTION','pairs':134,'joint_case_checks':joint_checks,'source_group_checks':checked_groups,'pair_mask_lineage_checks':source_masks,'extremum_provenance_checks':extremum_checks,'median_float_reconstruction_max_abs':max(median_float_differences),'summaries':summaries,'named_pair':named,'synthetic_checks':tests,'input_hashes_before_after':pins,'reviewer_sha256':sha(Path(__file__)),'production_imports':False,'mismatches':[],'scope':'Exposed-source exploratory stored-proxy result; not blind review, biological replication, calibrated error or normalization independence.'}
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'joint_cases':joint_checks,'medians_max_abs':max(median_float_differences),'named_pair':named,'synthetic_checks':len(tests),'receipt_sha256':sha(out/'verification.json')},indent=2))
if __name__=='__main__':main()
