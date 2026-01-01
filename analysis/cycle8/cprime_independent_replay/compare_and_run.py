"""Portable independent replay. Actual-source execution requires root authorization."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited before any read/output')
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as Q

def hashfile(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def compare(reference,production,eligible):
    """Compare every numeric record; validate arbitrary valid order certificates."""
    def require(test,msg):
        if not test:raise ArithmeticError(msg)
    def q(v):return None if v is None else Q(str(v))
    count={'source_rows':0,'state_cells':0,'group_gamma_records':0,'order_gamma_records':0,'edge_origins':0,'witness_values':0}
    rr={r['source_row']:r for r in reference['source_rows']};pp={r['source_row']:r for r in production['source_rows']}
    require(len(rr)==len(reference['source_rows']) and len(pp)==len(production['source_rows']) and set(rr)==set(pp),'Lost or duplicate source row')
    for k,r in rr.items():
        p=pp[k]
        for name in ['state','cell','intended','marker','exon','usable','source_score_consistent']:require(r[name]==p[name],f'Row {k} {name}')
        for name in ['DNA','RNA','source_score','ratio','score_difference']:require(q(r[name])==q(p[name]),f'Row {k} {name}')
        require(bool(r['missing_reason'])==bool(p['missing_reason']),f'Row {k} missing reason')
        count['source_rows']+=1
    if reference['status']=='SOURCE_REPRODUCTION_FAILED':
        require(production['status']=='SOURCE_REPRODUCTION_FAILED','Failed source reproduction concealed')
        return count
    rs={(r['state'],r['cell']):r for r in reference['state_cells']};ps={(r['state'],r['cell']):r for r in production['state_cells']}
    require(len(ps)==len(production['state_cells']) and set(rs)==set(ps),'State/cell loss or duplication')
    for k,r in rs.items():
        p=ps[k]
        for name in ['usable_rows','complete']:require(r[name]==p[name],f'State {k} {name}')
        require(sorted(r['source_rows'])==sorted(p['source_rows']),f'State {k} source lineage')
        for name in ['lower','upper','median']:require(q(r[name])==q(p[name]),f'State {k} {name}')
        count['state_cells']+=1
    rsen={q(s['gamma']):s for s in reference['sensitivities']};psen={q(s['gamma']):s for s in production['sensitivities']}
    require(set(rsen)==set(psen) and len(psen)==4,'Gamma loss/duplication')
    for gamma,s in rsen.items():
        prod=psen[gamma]
        rg={(g['intended'],g['cell']):g for g in s['group_cell_results']};pg={(g['intended'],g['cell']):g for g in prod['group_cell_results']}
        require(set(rg)==set(pg) and len(pg)==len(prod['group_cell_results']),'Group loss/duplication')
        for k,g in rg.items():
            p=pg[k]
            for name in ['marker_count','complete','maximum_attained','common_response_feasible']:require(g[name]==p[name],f'Group {k} {gamma} {name}')
            require(sorted(g['markers'])==sorted(p['markers']),f'Group {k} markers')
            for name in ['minimum_fold','maximum_fold','point_fold']:require(q(g[name])==q(p[name]),f'Group {k} {gamma} {name}')
            require(g['thresholds']==p['thresholds'],f'Group {k} grid')
            count['group_gamma_records']+=1
        ro={r['exon']:r for r in s['exon_shared_order_results']};po={r['exon']:r for r in prod['exon_shared_order_results']}
        require(set(ro)==set(po) and len(po)==len(prod['exon_shared_order_results']),'Exon loss/duplication')
        for exon,r in ro.items():
            p=po[exon];require(r['status']==p['status'] and r['state_count']==p['state_count'],f'Order {exon}/{gamma}')
            def edge_map(obj):
                out={}
                for edge in obj['edges']:
                    k=edge['lower_state'],edge['higher_state']
                    require(k not in out,'Duplicate order edge')
                    out[k]=sorted((a['context'],q(a['upper']),q(a['lower'])) for a in edge['origins'])
                return out
            re,pe=edge_map(r),edge_map(p);require(re==pe,f'Order edge/provenance mismatch {exon}/{gamma}')
            count['edge_origins']+=sum(len(v) for v in re.values())
            nodes={s for s,m in eligible['states'].items() if int(m['exon_number'])==exon}
            if r['status']=='feasible':
                order=p['latent_order'];require(set(order)==nodes and len(order)==len(nodes),'Invalid latent permutation')
                idx={n:i for i,n in enumerate(order)}
                require(all(idx[a]<idx[b] for a,b in re),'Non-topological latent order')
                require(set(p['responses'])==set(eligible['cell_contexts']),'Lost response context')
                for cell,vals in p['responses'].items():
                    require(set(vals)==nodes,'Lost response state')
                    values=[q(vals[n]) for n in order]
                    require(all(v>0 for v in values) and all(a<=b for a,b in zip(values,values[1:])),'Nonpositive/nonmonotone response')
                    for n,v in zip(order,values):
                        a=rs[n,cell];lo=q(a['lower'])/gamma;hi=None if a['upper'] is None else q(a['upper'])*gamma
                        require(v>=lo and (hi is None or v<=hi),'Response outside source interval')
                        count['witness_values']+=1
            else:
                cycle=list(p['cycle'])
                if len(cycle)>1 and cycle[0]==cycle[-1]:cycle.pop()
                require(len(cycle)>=2 and len(set(cycle))==len(cycle) and set(cycle)<=nodes,'Invalid cycle vertices')
                cycle_keys=list(zip(cycle,cycle[1:]+cycle[:1]));require(all(k in re for k in cycle_keys),'Unsupported cycle')
                cert={'edges':p['cycle_edges']};ce=edge_map(cert)
                require(set(ce)==set(cycle_keys) and all(ce[k]==re[k] for k in ce),'Cycle provenance mismatch')
            count['order_gamma_records']+=1
        require(s['summary']==prod['summary'],f'Summary mismatch {gamma}')
    return count

def main():
    parser=argparse.ArgumentParser()
    for key in ['source','protocol','eligibility','result','result-sha256','authorization','output']:parser.add_argument('--'+key,required=True)
    a=parser.parse_args();base=Path(__file__).parent
    seal=json.loads((base/'seal.json').read_text())
    for name,wanted in seal['files'].items():
        if hashfile(base/name)!=wanted:raise ValueError('Independent verifier byte drift: '+name)
    auth=json.loads(Path(a.authorization).read_text())
    if auth.get('purpose')!='B1_INDEPENDENT_REPLAY' or auth.get('authorized') is not True:raise PermissionError('Explicit root replay authorization required')
    if hashfile(a.protocol)!=seal['protocol_sha256'] or hashfile(a.eligibility)!=seal['eligibility_sha256']:raise ValueError('Frozen protocol/eligibility changed')
    if hashfile(a.source)!=seal['workbook_sha256'] or hashfile(a.result)!=a.result_sha256:raise ValueError('Source or result hash mismatch')
    if auth.get('protocol_sha256')!=seal['protocol_sha256'] or auth.get('source_workbook_sha256')!=seal['workbook_sha256'] or auth.get('production_result_sha256')!=a.result_sha256:raise PermissionError('Authorization artifact mismatch')
    output=Path(a.output)
    if output.exists():raise FileExistsError('Refusing existing output')
    # Only now is source response decoding permitted.
    import reference as core
    eligible=json.loads(Path(a.eligibility).read_text())
    rows,roundtrip=core.read_selected(a.source,eligible)
    if len(rows)!=396:raise ValueError('Unexpected selected outcome-row denominator')
    reconstructed=core.reconstruct(rows,eligible)
    production=json.loads(Path(a.result).read_text())
    coverage=compare(core.ready(reconstructed),production,eligible)
    output.mkdir(parents=True)
    for name,obj in [('independent_reconstruction.json',reconstructed),('openpyxl_roundtrip.json',roundtrip)]:
        (output/name).write_text(json.dumps(core.ready(obj),indent=2)+'\n',encoding='utf-8')
    receipt={'status':'PASS_ALL_INDEPENDENT_RECORDS','coverage':coverage,'openpyxl_roundtrip_differences':len(roundtrip),
             'arithmetic_uses':'Exact raw numeric XML decimal strings; no openpyxl rounding determines orders',
             'protocol_sha256':seal['protocol_sha256'],'source_sha256':seal['workbook_sha256'],
             'production_result_sha256':a.result_sha256,'historical_git_reauthenticated':False,
             'limitations':'Descriptive proxy replay, no biological calibration or independent experimental validation'}
    (output/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':main()
