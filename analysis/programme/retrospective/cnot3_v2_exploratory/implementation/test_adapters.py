"""Synthetic guard-only tests. No grouping, fitting or biological outcome access."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import pipeline as p
import predict as model


def fixture():
    ids=[f'ENSG{i:011d}' for i in range(5480)]
    before={key:(str(i)*64) for i,key in enumerate(('rows','protocol','original_gate','source_audit'),1)}
    old={'schema':'CNOT3_principal_reproduction_gate_v1','status':'held_for_source_diagnosis',
         'cohort_gene_count':5480,'cohort_GeneIDs':ids,'cohort_rows_sha256':before['rows'],
         'principal_discrepancies':['stored_HL_ratio_inconsistent:'+x for x in ids[:49]]+['published_cohort_size_mismatch']}
    gate={'schema':'CNOT3_operational_source_gate_v2','status':'operational_cohort_verified',
          'original_frozen_gate_status':'held_for_source_diagnosis','historical_published_cohort_reproduced':False,
          'operational_source_gate_passed':True,'unresolved_source_discrepancies':['published_cohort_size_mismatch'],
          'cohort_gene_count':5480,'cohort_rows_sha256':before['rows'],'prediction_protocol_sha256':before['protocol'],
          'original_frozen_gate_sha256':before['original_gate'],'independent_source_audit_sha256':before['source_audit'],
          'operational_checks':{key:True for key in p.OPERATIONAL_CHECKS},
          'cohort_GeneIDs':ids,'original_ratio_discrepancy_GeneIDs':ids[:49]}
    return gate,old,before


class AdapterTests(unittest.TestCase):
    def test_new_valid_gate_and_exact_eleven_checks(self):
        gate,old,before=fixture()
        self.assertEqual(len(p.OPERATIONAL_CHECKS),11)
        self.assertIn('no_top1000_boundary_tie',p.OPERATIONAL_CHECKS)
        self.assertIn('median_log2_change_top1000_greater_rest',p.OPERATIONAL_CHECKS)
        self.assertIsNone(p.verify_operational_gate(gate,old,before))

    def test_old_gate_and_false_pass_claim_are_rejected(self):
        gate,old,before=fixture()
        with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(old,old,before)
        cases=[('schema','CNOT3_principal_reproduction_gate_v1'),('status','passed'),
               ('original_frozen_gate_status','passed'),('historical_published_cohort_reproduced',True),
               ('reproduction_gate_passed',True),('operational_source_gate_passed',False)]
        for key,value in cases:
            with self.subTest(key=key):
                changed=copy.deepcopy(gate);changed[key]=value
                with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(changed,old,before)

    def test_missing_false_unknown_and_old_checks_rejected(self):
        gate,old,before=fixture()
        for key in p.OPERATIONAL_CHECKS:
            changed=copy.deepcopy(gate);changed['operational_checks'][key]=False
            with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(changed,old,before)
            changed=copy.deepcopy(gate);del changed['operational_checks'][key]
            with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(changed,old,before)
        changed=copy.deepcopy(gate);changed['operational_checks']['count_5481']=True
        with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(changed,old,before)

    def test_hashes_count_unresolved_and_49_retention(self):
        gate,old,before=fixture()
        for key in ('cohort_rows_sha256','prediction_protocol_sha256','original_frozen_gate_sha256','independent_source_audit_sha256'):
            changed=copy.deepcopy(gate);changed[key]='f'*64
            with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(changed,old,before)
        for key,value in [('cohort_gene_count',5481),('unresolved_source_discrepancies',[]),
                          ('unresolved_source_discrepancies',['published_cohort_size_mismatch','new_unexplained_failure']),
                          ('original_ratio_discrepancy_GeneIDs',gate['original_ratio_discrepancy_GeneIDs'][:-1]),
                          ('original_ratio_discrepancy_GeneIDs',gate['original_ratio_discrepancy_GeneIDs'][1:]+[gate['cohort_GeneIDs'][50]])]:
            changed=copy.deepcopy(gate);changed[key]=value
            with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(changed,old,before)
        changed=copy.deepcopy(gate);changed['cohort_GeneIDs'][0]='ENSG99999999999'
        with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(changed,old,before)

    def test_original_gate_cannot_be_substituted(self):
        gate,old,before=fixture()
        for key,value in [('status','passed'),('cohort_gene_count',5481),('cohort_rows_sha256','e'*64),
                          ('principal_discrepancies',old['principal_discrepancies']+['other_error'])]:
            changed=copy.deepcopy(old);changed[key]=value
            with self.assertRaises(p.c.IntegrityError):p.verify_operational_gate(gate,changed,before)

    def test_source_audit_bytes_pinned_before_any_grouping(self):
        # Only fake metadata files. The mutated audit hash must fail before
        # reference pins, references, grouping or model data can be processed.
        with tempfile.TemporaryDirectory() as td:
            root=Path(td);gate,old,_=fixture()
            paths={name:root/(name+'.json') for name in ('rows','original_ids','features','metadata','fasta','upstream_gate','protocol','original_gate','source_audit')}
            for path in paths.values():path.write_text('{}',encoding='utf-8')
            rows=[dict(GeneID=g,HL_WT=1.,HL_CNOT3KO=2.,LOG2_FC_IP=0.,SE_IP=1.) for g in gate['cohort_GeneIDs']]
            paths['rows'].write_text(json.dumps(rows),encoding='utf-8')
            old['cohort_rows_sha256']=p.sha(paths['rows'])
            paths['original_gate'].write_text(json.dumps(old),encoding='utf-8')
            for key,role in [('cohort_rows_sha256','rows'),('prediction_protocol_sha256','protocol'),('original_frozen_gate_sha256','original_gate'),('independent_source_audit_sha256','source_audit')]:gate[key]=p.sha(paths[role])
            paths['upstream_gate'].write_text(json.dumps(gate),encoding='utf-8')
            auth={'operation':'CNOT3_operational_cohort_group_partition_v2','real_execution_authorized':True,
                  'operational_source_gate_passed':True,'input_sha256':{k:p.sha(v) for k,v in paths.items()},
                  'code_sha256':{'cohort.py':p.sha(p.HERE/'cohort.py'),'pipeline.py':p.sha(p.__file__)}}
            authorization=root/'authorization.json';authorization.write_text(json.dumps(auth),encoding='utf-8')
            paths['source_audit'].write_text('{"changed":true}',encoding='utf-8')
            with mock.patch.object(p.c,'group_proteins',side_effect=AssertionError('Grouping must not run')),mock.patch.object(p,'references',side_effect=AssertionError('References must not open')):
                with self.assertRaisesRegex(p.c.IntegrityError,'hashes differ'):
                    p.execute(paths,authorization,root/'out')
            self.assertFalse((root/'out').exists())

    def test_prediction_guard_accepts_only_explicit_exploratory_semantics(self):
        # These are synthetic identity/group metadata, with no target or feature data.
        rows=[dict(GeneID=f'SYN{i:05d}',group_ID=f'SYN{i:05d}') for i in range(5480)]
        train=[r for r in rows if model.group_fold(r['group_ID'],'outer')!=0]
        holdout=[r for r in rows if model.group_fold(r['group_ID'],'outer')==0]
        receipt=dict(status='root_authorized_operational_prediction',analysis_status='post_exposure_exploratory',
                     historical_published_cohort_reproduced=False,original_frozen_gate_status='held_for_source_diagnosis',
                     operational_source_gate_passed=True,reference_eligibility_passed=True,grouping_complete=True,
                     operational_source_gate_sha256='1'*64,original_frozen_gate_sha256='2'*64,independent_source_audit_sha256='3'*64,
                     operational_slam_overlap_count=5480,eligible_count=5480)
        got=model.check_design(train,holdout,receipt)
        self.assertEqual(got['status'],'root_authorized_operational_prediction')
        for key,value in [('status','root_authorized_prediction'),('reproduction_gate_passed',True),
                          ('historical_published_cohort_reproduced',True),('analysis_status','validation'),
                          ('operational_source_gate_passed',False),('operational_source_gate_sha256','unverified')]:
            changed=copy.deepcopy(receipt);changed[key]=value
            with self.assertRaises(ValueError):model.check_design(train,holdout,changed)

    def test_original_hash_salts_and_thresholds_preserved(self):
        import hashlib
        for name in ('outer','inner'):
            expected=int.from_bytes(hashlib.sha256(('cnot3-focus-v1|'+name+'|SYN00001').encode()).digest(),'big')%5
            self.assertEqual(model.group_fold('SYN00001',name),expected)
        self.assertEqual(model.MINIMUMS,dict(training_genes=2000,holdout_genes=300,training_groups=100,holdout_groups=30,assessment_genes=100,inner_training_groups=50))
        self.assertEqual(p.c.LIMITS.posting_pair_events,20_000_000)
        self.assertEqual(p.c.LIMITS.candidate_pairs,50_000)
        self.assertEqual(p.c.LIMITS.alignment_cells,3_000_000_000)
        self.assertEqual(p.c.LIMITS.elapsed_seconds,300.)


if __name__=='__main__':unittest.main(verbosity=2)
