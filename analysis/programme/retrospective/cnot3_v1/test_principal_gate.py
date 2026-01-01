"""Synthetic integrity and primary-versus-secondary gate counterexamples."""
import copy
import importlib.util
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('gate',HERE/'principal_gate.py')
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)

class GateTests(unittest.TestCase):
    def setUp(self):
        ids=['ENSG'+str(i).zfill(11) for i in range(5481)]
        self.rows=[dict(GeneID=i,HL_WT=1.0,HL_CNOT3KO=2.0,LOG2_FC_IP=0.0,SE_IP=1.0) for i in ids]
        self.r={'n_cohort':5481,'n_high':1000,'n_remainder':4481,'cohort_ids':ids,'discrepancies':[],
          'half_life_ratio_diagnostics':[{'GENE_ID':i,'consistent':True} for i in ids],
          'IP_statistic_diagnostics':[{'GENE_ID':i,'consistent':True} for i in ids],
          'selection_ties':{'status':'selected','boundary_tie_crosses':False},
          'group_effects':{'HL_FC_LOG2':{'median_difference':.1}},'unresolved_source_rules':['Synthetic only']}
    def test_declared_secondary_warnings_do_not_change_primary_gate(self):
        self.assertEqual(g.assess(self.r,self.rows)['status'],'passed')
        for warning in g.SECONDARY:
            x=copy.deepcopy(self.r);x['discrepancies']=[warning]
            self.assertEqual(g.assess(x,self.rows)['status'],'passed')
    def test_unknown_and_nonpositive_precision_issues_stop(self):
        for issue in ['unknown_issue','negative_SE_IP:outside_final_cohort','zero_SE_IP:ENSG00000000000']:
            x=copy.deepcopy(self.r);x['discrepancies']=[issue]
            self.assertEqual(g.assess(x,self.rows)['status'],'held_for_source_diagnosis')
    def test_arithmetic_tie_direction_and_membership_failures(self):
        variants=[]
        x=copy.deepcopy(self.r);x['half_life_ratio_diagnostics'][0]['consistent']=False;variants.append(x)
        x=copy.deepcopy(self.r);x['selection_ties']['boundary_tie_crosses']=True;variants.append(x)
        x=copy.deepcopy(self.r);x['group_effects']['HL_FC_LOG2']['median_difference']=0;variants.append(x)
        x=copy.deepcopy(self.r);x['cohort_ids']=x['cohort_ids'][1:];variants.append(x)
        for x in variants:self.assertEqual(g.assess(x,self.rows)['status'],'held_for_source_diagnosis')

if __name__=='__main__':unittest.main(verbosity=2)
