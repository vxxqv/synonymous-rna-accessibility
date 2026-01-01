"""Verify corrected robustness revision using the held independent weak-order oracle."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import ast
import difflib
import hashlib
import json
import subprocess
import sys
import time
import review_synthetic as reference

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
EXPECTED='7a0798391e1987f06dbb2d09f6771aa10ffa39431d2f765400e47e594a9f395e'
PROTO='6d3edb6f28f8a71a5cf951c7163e4775ceff371c74071633d82fbe32e967dc58'
need=reference.need
sha=reference.sha

def scrub(obj):
    if isinstance(obj,dict):return {k:scrub(v) for k,v in obj.items() if not k.startswith('decimal_')}
    if isinstance(obj,list):return [scrub(v) for v in obj]
    return obj

def main():
    dest=BASE/'delta_receipt.json';need(not dest.exists(),'Refusing existing receipt');start=time.perf_counter()
    source=ROOT/'work/research_cycle7/robustness.py';raw=source.read_bytes()
    proto=ROOT/'work/research_cycle7/robustness_protocol.md';proto_raw=proto.read_bytes()
    need(sha(raw)==EXPECTED and sha(proto_raw)==PROTO,'Candidate drift')
    snapshot=BASE/'corrected_robustness.py';need(not snapshot.exists(),'Snapshot exists');snapshot.write_bytes(raw)
    (BASE/'corrected_protocol.md').write_bytes(proto_raw)
    original=(BASE/'original_robustness.py').read_bytes()
    (BASE/'correction.diff').write_text(''.join(difflib.unified_diff(original.decode().splitlines(True),raw.decode().splitlines(True),fromfile='original_robustness.py',tofile='corrected_robustness.py')),encoding='utf-8')
    author=reference.load(snapshot,'corrected_robustness');old=reference.load(BASE/'original_robustness.py','original_robustness')
    core=reference.load(reference.CORE,'frozen_core_delta')
    X=[[1,0,0],[1,1,0],[1,0,1],[1,1,1]];orders=reference.weak_orders()
    counts={'finite':0,'infinite':0,'point_systems':0,'interval_systems':0}
    for vals in product((Q(0),Q(1,2),Q(1)),repeat=4):
        for mask in product((0,1),repeat=4):
            for metric in ('additive','multiplicative_squared'):
                status=reference.audit_one(author,core,X,orders,list(vals),list(vals),mask,metric)
                counts[status]+=1;counts['point_systems']+=1
    alphabet=[(Q(0),Q(0)),(Q(1,2),Q(1)),(Q(0),None)]
    for intervals in product(alphabet,repeat=4):
        L=[a for a,b in intervals];U=[b for a,b in intervals]
        for mask in ((0,1,1,1),(0,0,0,0),(1,1,1,1)):
            for metric in ('additive','multiplicative_squared'):
                status=reference.audit_one(author,core,X,orders,L,U,mask,metric)
                counts[status]+=1;counts['interval_systems']+=1
    # These masks must fail before an optimizer can run.
    class Trap:
        def order_certificate(self,*args):raise AssertionError('Core called before input guard')
        def features(self,*args):raise AssertionError('Features called before input guard')
    for mask in ([0.0,1.0],[False,True],[0,2]):
        try:author.radius(Trap(),[[1],[1]],[Q(1),Q(3,5)],[Q(1),Q(3,5)],mask,'multiplicative_squared')
        except ValueError as error:need('Integer binary' in str(error),'Wrong mask failure')
        else:raise AssertionError('Unsafe mask accepted')
    exact=author.radius(core,[[1],[1]],[Q(1),Q(3,5)],[Q(1),Q(3,5)],[0,1],'multiplicative_squared')
    need(exact['parameter_exact']=='25/9','Exact integer-mask counterpart')
    states=list(product(range(3),range(3),range(2),range(2),range(2)))
    data={'vertices_data':[{'vertex':i,'codon_states':list(s),'lower':'1','upper':'1','scenario_ratios':['1']*27} for i,s in enumerate(states)]}
    actual=author.analyze(core,data);previous=old.analyze(core,data)
    need(scrub(actual)==scrub(previous),'Valid deletion/scientific payload changed')
    need(len(actual['all_single_vertex_deletions'])==72 and len(actual['all_single_preparation_deletions'])==12,'Deletion coverage')
    for omit,row in enumerate(actual['all_single_vertex_deletions']):
        need(row['omitted_vertex']==omit and row['retained_vertex_ids']==[i for i in range(72) if i!=omit],'Row deletion identities')
    for index,row in enumerate(actual['all_single_preparation_deletions']):
        kind='output' if index<9 else 'input';omitted=index if index<9 else index-9
        columns=[3*b+t for b in range(9) for t in range(3) if (b if kind=='output' else t)!=omitted]
        need(row['kind']==kind and row['omitted_zero_based']==omitted and row['retained_scenarios']==columns,'Preparation identities')
    input_guards=[]
    for name in ('unknown_upper','unknown_scenario','short_scenario','permuted_vertex'):
        bad=json.loads(json.dumps(data))
        if name=='unknown_upper':bad['vertices_data'][0]['upper']=None
        elif name=='unknown_scenario':bad['vertices_data'][0]['scenario_ratios'][0]=None
        elif name=='short_scenario':bad['vertices_data'][0]['scenario_ratios'].pop()
        else:bad['vertices_data'][0]['vertex']=1
        try:author.analyze(Trap(),bad)
        except ValueError:input_guards.append(name)
        else:raise AssertionError('Missing domain guard')
    root_cases=(Q(0),Q(1),Q(2),Q(4),Q(25,9),Q(10**40+1))
    for value in root_cases:
        bounds=author.sqrt_enclosure(value)
        lo,hi=Q(bounds['lower']),Q(bounds['upper'])
        need(lo*lo<=value<=hi*hi and Q(bounds['decimal_lower_exact'])==lo and Q(bounds['decimal_upper_exact'])==hi,'Exact decimal enclosure')
    # Current guard tests use only {} as a fake result, never a biological file.
    marker=BASE/'corrected_core_executed.txt';poison=BASE/'corrected_synthetic_core.py'
    poison.write_text('from pathlib import Path\nPath('+repr(str(marker))+').write_text("executed")\n',encoding='utf-8')
    toy=BASE/'toy_result.json';toy.write_bytes(b'{}')
    freeze=BASE/'toy_freeze.json'
    frozen={'files':{'corrected_robustness.py':EXPECTED},'primary_result_sha256':sha(b'{}'),'ordinal_core_sha256':'0'*64}
    freeze.write_text(json.dumps(frozen),encoding='utf-8')
    guard_checks=[]
    def invoke(label,extra=(),freeze_path=freeze,out=None):
        output=BASE/(label+'.json') if out is None else out
        p=subprocess.run([sys.executable,*extra,str(snapshot),'--core',str(poison),'--out',str(output),
                          '--result',str(toy),'--freeze',str(freeze_path)],capture_output=True,text=True)
        need(p.returncode!=0 and not marker.exists(),'Unverified core executed: '+label)
        if out is None:need(not output.exists(),'Failure wrote output')
        guard_checks.append({'case':label,'returncode':p.returncode,'last_error':p.stderr.splitlines()[-1]})
    invoke('reject_wrong_core_hash')
    invoke('reject_missing_freeze',freeze_path=BASE/'nonexistent_freeze_v2.json')
    frozen['primary_result_sha256']='0'*64;freeze.write_text(json.dumps(frozen),encoding='utf-8')
    invoke('reject_wrong_result_hash')
    invoke('reject_optimization_O',extra=('-O',));invoke('reject_optimization_OO',extra=('-OO',))
    occupied=BASE/'existing_output.json';occupied.write_text('original',encoding='utf-8')
    invoke('reject_overwrite',out=occupied);need(occupied.read_text()=='original','Overwrite changed content')
    oldfuncs={n.name:ast.dump(n,include_attributes=False) for n in ast.parse(original.decode()).body if isinstance(n,ast.FunctionDef)}
    newfuncs={n.name:ast.dump(n,include_attributes=False) for n in ast.parse(raw.decode()).body if isinstance(n,ast.FunctionDef)}
    unchanged=[n for n in oldfuncs if oldfuncs[n]==newfuncs[n]]
    need('fields' in unchanged and 'breakpoints' in unchanged,'Mathematical transformation changed')
    need(source.read_bytes()==raw and proto.read_bytes()==proto_raw,'Source drift during audit')
    result={'status':'PASS_FOR_DECLARED_COMPLETE_OBSERVED_DOMAIN', 'source_sha256':EXPECTED,'protocol_sha256':PROTO,
            'original_sha256':sha(original),'core_sha256':sha(reference.CORE.read_bytes()),
            'independent_weak_orders':17,'independent_counts':counts,
            'all72_row_and12_preparation_deletions_verified':True,
            'valid_synthetic_deletion_payload_equal_except_decimal_display_fields':True,
            'mask_rejections':3,'exact_counterpart_squared_radius':'25/9','domain_guards':input_guards,
            'exact_decimal_root_cases':6,'preimport_and_output_guards':guard_checks,
            'unchanged_function_asts':unchanged,'real_outcomes_read':False,
            'runtime_seconds':time.perf_counter()-start,'verifier_sha256':sha(Path(__file__).read_bytes()),
            'limitations':['Generic radius supports unknown upper intervals; observed adapter intentionally rejects incomplete fields before computation.',
                           'Floating LP proposes certificates; exact primal/dual checks remain mandatory and solver failures must stop.',
                           'Provided freeze is a trusted external input, not proof by itself of human authorization or Git publication.',
                           'The finite Q2 oracle validates this test domain; general correctness relies on the reviewed breakpoint proof.']}
    dest.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','source_sha256','independent_counts','runtime_seconds','real_outcomes_read')}))

if __name__=='__main__':main()
