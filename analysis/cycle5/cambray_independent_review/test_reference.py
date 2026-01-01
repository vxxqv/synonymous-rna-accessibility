"""Synthetic-only independent tests. Never open the deposited phenotype source."""
import sys
if not __debug__:raise RuntimeError('Optimization refused before input access')
sys.dont_write_bytecode=True
import copy
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import unittest
import reference as r

class Tests(unittest.TestCase):
    def test_all16_patterns(self):
        positives=0
        for a,b,c,d in product((False,True),repeat=4):
            vals={0:(1,1),1:(int(a),int(c)),2:(int(b),int(d)),3:(1,1)}
            result=r.classify(vals,r.SQUARE_EDGES,0,3,(1,1))
            expected=(a or b) and(c or d) and not((a and c) or(b and d))
            self.assertEqual(result['status']=='joint_only',expected)
            positives+=expected
        self.assertEqual(positives,2)

    def test_additive_antagonism_aligned(self):
        a=[1,0,2,1];b=[1,2,0,1]
        vals={i:(a[i],b[i]) for i in range(4)}
        self.assertEqual(r.classify(vals,r.SQUARE_EDGES,0,3)['status'],'joint_only')
        self.assertEqual(r.interaction(a)['contrast'],0)
        self.assertEqual(r.interaction(b)['contrast'],0)
        self.assertEqual(r.classify({i:(a[i],a[i]) for i in range(4)},r.SQUARE_EDGES,0,3)['status'],'joint_connected')

    def test_ties_and_unacceptable_endpoint(self):
        vals={0:(2,2),1:(1,0),2:(0,1),3:(2,2)}
        self.assertEqual(r.classify(vals,r.SQUARE_EDGES,0,3,(1,1))['status'],'joint_only')
        self.assertEqual(r.classify(vals,r.SQUARE_EDGES,0,3,(3,1))['status'],'ineligible_endpoint')
        self.assertEqual(r.classify(vals,r.SQUARE_EDGES,0,3,(1,1),(1e-10,0))['status'],'joint_only')

    def test_nonmonotone_detour(self):
        vertices=list(product((0,1),repeat=3));edges=[(a,b) for a in vertices for b in vertices if a<b and sum(x!=y for x,y in zip(a,b))==1]
        vals={x:(x[0]-x[1]+x[2],-x[0]+x[1]+x[2]) for x in vertices}
        self.assertEqual(r.classify(vals,edges,(0,0,0),(1,1,0))['status'],'joint_connected')
        local={x:v for x,v in vals.items() if x[2]==0};local_edges=[(a,b) for a,b in edges if a in local and b in local]
        self.assertEqual(r.classify(local,local_edges,(0,0,0),(1,1,0))['status'],'joint_only')

    def test_interval_certificate_and_nonnesting(self):
        ints={0:((0,0),(0,0)),1:((1,2),(-2,-1)),2:((-2,-1),(1,2)),3:((0,0),(0,0))}
        self.assertTrue(r.bounds(ints,r.SQUARE_EDGES,0,3,(0,0))['joint_only_certificate'])
        reps={0:[0,10,0,10,0,10,0],1:[0]*4+[-1]*3,2:[-1]*4+[0]*3,3:[0,10,0,10,0,10,0]}
        out=r.preparation_summary(reps)
        self.assertTrue(out['range_certificate']['joint_only_certificate'])
        self.assertFalse(out['all12_joint_only'])
        self.assertEqual(len(out['preparation_statuses']),12)

    def test_missing_trait_known_failure(self):
        vs=[(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(1,1,1)]
        edges=[(a,b) for a in vs for b in vs if a<b and sum(x!=y for x,y in zip(a,b))==1]
        ints={x:((0,0),(0,0)) for x in vs}
        ints[(1,0,0)]=((1,1),(-1,-1));ints[(0,1,0)]=((-1,-1),(1,1))
        ints[(1,0,1)]=(None,(0,0))
        result=r.bounds(ints,edges,(0,0,0),(1,1,0),(0,0))
        self.assertTrue(result['joint_optimistic']);self.assertFalse(result['joint_pessimistic'])
        ints[(1,0,1)]=(None,(-1,-1))
        self.assertTrue(r.bounds(ints,edges,(0,0,0),(1,1,0),(0,0))['joint_only_certificate'])

    def test_missing_is_not_zero_and_preparation_structure(self):
        reps={i:[1.]*7 for i in range(4)};reps[1][0]=None
        self.assertEqual(r.preparation_summary(reps)['status'],'incomplete_required_measurements')
        vals={0:(1,1),1:(None,1),2:(0,1),3:(1,1)}
        self.assertEqual(r.classify(vals,r.SQUARE_EDGES,0,3)['status'],'incomplete_point_field')

    def test_local_prediction_exact_and_key_guard(self):
        a=[2,5,7,10];b=[5,3,9,7]
        for target in range(4):
            known={i:(a[i],b[i]) for i in range(4) if i!=target}
            self.assertEqual(r.local_predictions(known,target)['additive'],[a[target],b[target]])
        class Poison:
            def __len__(self):raise AssertionError('Hidden value touched')
        with self.assertRaises(ValueError):r.local_predictions({0:Poison(),1:Poison(),2:Poison(),3:Poison()},3)
        known={0:(5,6),1:(0,0),2:(0,0)}
        self.assertEqual(r.local_predictions(known,3)['additive'],[-5,-6])

    def test_interaction_projector(self):
        y=[1,4,3,10];out=r.interaction(y)
        self.assertEqual(out['contrast'],4)
        self.assertEqual(out['projection'],[0,5,4,9])
        self.assertAlmostEqual(r.interaction(out['projection'])['contrast'],0)
        self.assertIsNotNone(out['natural_log'])
        self.assertIsNone(r.interaction([0,1,1,2])['natural_log'])
        # Three observations do not determine the fourth interaction contrast.
        self.assertNotEqual(r.interaction([1,2,3,4])['contrast'],r.interaction([1,2,3,5])['contrast'])

    def test_critical_cells_exhaustive_small_values(self):
        count=0
        for a,b,c,d in product((-1.,0.,1.),repeat=4):
            vals={0:(1.,1.),1:(a,c),2:(b,d),3:(1.,1.)}
            regions=r.critical_regions(vals,r.SQUARE_EDGES,0,3)
            for x,y in product((-3.,-1.,-.5,0.,.5,1.),repeat=2):
                matches=[cell for cell in regions if all((lo is None or z>lo) and z<=hi for z,(lo,hi) in zip((x,y),cell['cells']))]
                self.assertEqual(len(matches),1)
                self.assertEqual(matches[0]['classification'],r.classify(vals,r.SQUARE_EDGES,0,3,(x,y))['status'])
                count+=1
        self.assertEqual(count,2916)

    def make_cases(self,n=10):
        return [{'group':str(i%8),'truth':[1.,2.],
                 'predictions':{'additive':[1.1,2.1],'neighbor_mean':[1.3,2.3],'three_mean':[1.5,2.5]},
                 'true_class':bool(i%2),'predicted_class':{m:bool(i%2) for m in r.METHODS}} for i in range(n)]

    def test_group_metrics_success_failures(self):
        cases=self.make_cases();out=r.grouped_metrics(cases)
        self.assertTrue(out['success']);self.assertEqual(out['n_groups'],8)
        self.assertAlmostEqual(out['metrics']['additive']['MAE'][0],.1)
        self.assertAlmostEqual(out['metrics']['additive']['RMSE'][1],.1)
        cases[0]['predicted_class']['additive']=not cases[0]['true_class']
        self.assertFalse(r.grouped_metrics(cases)['success'])
        one=self.make_cases()
        for row in one:row['true_class']=False
        self.assertFalse(r.grouped_metrics(one)['success'])
        zero=self.make_cases()
        for row in zero:row['predictions']['neighbor_mean']=row['truth'][:]
        self.assertFalse(r.grouped_metrics(zero)['success'])
        self.assertFalse(r.grouped_metrics(self.make_cases(9))['coverage_pass'])
        self.assertFalse(r.grouped_metrics([])['success'])

    def test_group_balance_not_case_balance(self):
        cases=self.make_cases(3)
        for i,row in enumerate(cases):
            row['group']='a' if i<2 else 'b'
            row['predictions']['additive']=[1.,2.] if i<2 else [3.,4.]
        out=r.grouped_metrics(cases)['metrics']['additive']
        self.assertEqual(out['MAE'],[1.,1.]);self.assertEqual(out['RMSE'],[math.sqrt(2),math.sqrt(2)])

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Tests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)
