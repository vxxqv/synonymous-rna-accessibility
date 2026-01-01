"""Read authorized saved cPRIME result; summarize scope, not numerical replay."""
if not __debug__:
    raise RuntimeError('Optimized execution is forbidden.')

from pathlib import Path
from fractions import Fraction as F
import hashlib
import json

HERE = Path(__file__).resolve().parent
CYCLE = HERE.parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    result_path=CYCLE/'cprime_results/primary_v1/result.json'
    assert sha(result_path)=='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf'
    result=json.loads(result_path.read_text())
    freeze_path=CYCLE/'cprime_pilot/freeze.json'
    freeze=json.loads(freeze_path.read_text())
    assert result['freeze_sha256']==sha(freeze_path)
    assert result['status']=='COMPLETE_DESCRIPTIVE_AUDIT'
    assert len(result['source_rows'])==396 and len(result['state_cells'])==132
    eligibility=json.loads((CYCLE/'cprime_pilot/eligibility.json').read_text())
    for name,pin in freeze['files'].items():assert sha(CYCLE/'cprime_pilot'/name)==pin
    assert len(eligibility['groups'])==30 and len(eligibility['states'])==66
    assert all(x['complete'] for x in result['state_cells'])
    observations={(x['state'],x['cell']):x for x in result['state_cells']}
    extracted=[]
    for s in result['sensitivities']:
        assert len(s['group_cell_results'])==60 and len(s['exon_shared_order_results'])==5
        assert sorted(x['exon'] for x in s['exon_shared_order_results'])==[17,25,27,39,42]
        cycles=[]
        for x in s['exon_shared_order_results']:
            if x['status']!='infeasible':continue
            assert s['gamma']=='1'
            states=x['cycle']
            assert states[0]==states[-1]
            assert all(int(eligibility['states'][v]['exon_number'])==x['exon'] for v in states)
            for e in x['cycle_edges']:
                assert e['origins']
                for origin in e['origins']:
                    a=observations[e['lower_state'],origin['context']]
                    b=observations[e['higher_state'],origin['context']]
                    assert F(origin['upper'])==F(a['upper'])<F(b['lower'])==F(origin['lower'])
            intended=sorted({eligibility['states'][v]['var_id'] for v in states})
            cycles.append({'exon':x['exon'],'cycle_edges':len(states)-1,'distinct_intended_labels':len(intended),
                           'within_one_intended_label':len(intended)==1,'intended_labels':intended,
                           'certificate_origins_match_primary_state_intervals':True})
        extracted.append({'gamma':s['gamma'],'reported_summary':s['summary'],
                          'reported_fixed_grid':{t:{'lower_at_least':sum(g['thresholds'][t]['minimum_at_least'] for g in s['group_cell_results']),
                                                   'upper_at_most':sum(g['thresholds'][t]['maximum_at_most'] for g in s['group_cell_results'])} for t in ['5/4','3/2','2','4']},
                          'cycle_scope':cycles})
    files=['cprime_pilot/protocol.md','cprime_pilot/freeze.json','cprime_pilot/eligibility.json','cprime_results/primary_v1/result.json','cprime_results/primary_v1/receipt.json','cprime_gate/report.md','cprime_gate/measurement_contract.json','cprime_gate/marker_context_primary.txt','cprime_gate/article.xml']
    packet={'status':'interpretation_evidence_extracted','input_sha256':{n:sha(CYCLE/n) for n in files},
            'public_freeze_commit_parent_reported':'b6f83e550a8e6f7bbf7901b9e9d69e01956eaac5',
            'historical_remote_order_independently_verified_here':False,'source_counts':{'intended_labels':30,'marker_states':66,'state_cells':132,'group_cells':60,'source_rows':396,'complete_state_cells':132},
            'sensitivities':extracted,'scope':'Saved summaries and exact cycle-support/scope checks only. A separately assigned agent performs full independent numerical replay. No raw workbook, new endpoint or COMT outcome access.'}
    (HERE/'evidence.json').write_text(json.dumps(packet,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':packet['status'],'source_counts':packet['source_counts'],'sensitivities':extracted},indent=2))


if __name__=='__main__':main()
