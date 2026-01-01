"""Write exact diff and AST evidence, then run synthetic guard tests only."""
import ast
import difflib
import hashlib
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parent
ALLOWED={'cohort.py':{'design_guards'},'pipeline.py':{'execute','main'},
         'predict.py':{'check_design','run_predict','run_evaluate'}}
ADDED={'cohort.py':set(),'pipeline.py':{'verify_operational_gate'},'predict.py':set()}
CONST_CHANGED={'cohort.py':set(),'pipeline.py':{'PRINCIPAL_CHECKS','SECONDARY_WARNINGS','OPERATIONAL_CHECKS'},'predict.py':{'VERSION'}}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def node_hash(node):return hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest()


def audit():
    receipts={};diff=[]
    for filename in ALLOWED:
        original=ROOT/'originals'/filename;candidate=ROOT/filename
        old=original.read_text(encoding='utf-8');new=candidate.read_text(encoding='utf-8')
        oldtree,newtree=ast.parse(old),ast.parse(new)
        functions=lambda tree:{n.name:n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef,ast.ClassDef))}
        of,nf=functions(oldtree),functions(newtree)
        changed={k for k in of.keys()&nf.keys() if node_hash(of[k])!=node_hash(nf[k])}
        if changed!=ALLOWED[filename] or set(nf)-set(of)!=ADDED[filename] or set(of)-set(nf):
            raise AssertionError(('Unexpected AST function change',filename,changed))
        assignments=lambda tree:{n.targets[0].id:n for n in tree.body if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)}
        oa,na=assignments(oldtree),assignments(newtree)
        constant_changes={k for k in oa.keys()|na.keys() if k not in oa or k not in na or node_hash(oa[k])!=node_hash(na[k])}
        if constant_changes!=CONST_CHANGED[filename]:raise AssertionError(('Unexpected constant change',filename,constant_changes))
        unchanged={k:node_hash(of[k]) for k in sorted(of) if k not in changed}
        receipts[filename]={'original_sha256':sha(original),'candidate_sha256':sha(candidate),
                            'identical_function_or_class_AST':unchanged,'changed_guard_functions':sorted(changed),
                            'new_guard_functions':sorted(ADDED[filename]),'changed_label_or_gate_constants':sorted(constant_changes)}
        diff.extend(difflib.unified_diff(old.splitlines(keepends=True),new.splitlines(keepends=True),
                                         fromfile='original_v1/'+filename,tofile='operational_exploratory_v2/'+filename))
    (ROOT/'exact_changes.diff').write_text(''.join(diff),encoding='utf-8',newline='\n')
    (ROOT/'ast_evidence.json').write_text(json.dumps(receipts,indent=2)+'\n',encoding='utf-8')
    suite=unittest.defaultTestLoader.discover(str(ROOT),pattern='test_adapters.py')
    with (ROOT/'test_log.txt').open('w',encoding='utf-8') as handle:
        result=unittest.TextTestRunner(stream=handle,verbosity=2).run(suite)
    report={'status':'ready_for_independent_review' if result.wasSuccessful() else 'failed_guard_tests',
            'tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
            'actual_grouping_runs':0,'actual_model_fits':0,'outcome_reads':0,
            'function_change_contract':ALLOWED,'source_sha256':{p.name:sha(p) for p in ROOT.glob('*.py')},
            'artifact_sha256':{name:sha(ROOT/name) for name in ('exact_changes.diff','ast_evidence.json','test_log.txt')},
            'protocol_sha256':sha(ROOT.parent/'protocol.md')}
    report['function_change_contract']={key:sorted(value) for key,value in ALLOWED.items()}
    (ROOT/'preparation_receipt.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':audit()
