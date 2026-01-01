"""Disclosed post-access schema adapter; sealed numerical code stays unchanged.

The initial actual comparison stopped because its independent summary stored an
infeasible-exon count while production stores exon IDs. Derive those IDs from
the independently reconstructed order records. No numeric formula is changed.
"""
if not __debug__:raise RuntimeError('Optimized Python prohibited')
from pathlib import Path
import hashlib,json,sys
import compare_and_run as original

base=Path(__file__).parent
original_compare=original.compare

def adapted_compare(reference,production,eligible):
    if reference['status']!='SOURCE_REPRODUCTION_FAILED':
        for sensitivity in reference['sensitivities']:
            ids=sorted(r['exon'] for r in sensitivity['exon_shared_order_results'] if r['status']=='infeasible')
            if sensitivity['summary']['shared_order_infeasible_exons']!=len(ids):raise ArithmeticError('Independent count/list inconsistency')
            sensitivity['summary']['shared_order_infeasible_exons']=ids
    return original_compare(reference,production,eligible)

original.compare=adapted_compare
if __name__=='__main__':
    original.main()
    output=Path(sys.argv[sys.argv.index('--output')+1])
    p=output/'verification.json';d=json.loads(p.read_text())
    d['schema_adapter']={'path':Path(__file__).name,'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                         'change':'Only independent summary count-to-exon-ID-list shape for comparison; original numerical core and reconstruction remain unchanged',
                         'initial_failure_log_sha256':hashlib.sha256((base/'actual_v1_execution.log').read_bytes()).hexdigest()}
    p.write_text(json.dumps(d,indent=2)+'\n')
