"""Synthetic-only verification with independent exhaustive simple-path oracles."""
import hashlib
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path
import unittest
import numpy as np

import graph_rules as rules

SQUARE = [[1,2],[0,3],[0,3],[1,2]]
CUBE = [[v^1,v^2,v^4] for v in range(8)]
ROOT = Path(__file__).resolve().parent


def simple_paths(graph, start, target):
    """Direct path enumeration, independent from production graph search."""
    def visit(path):
        if path[-1] == target:
            yield path
            return
        for vertex in graph[path[-1]]:
            if vertex not in path:
                yield from visit(path+[vertex])
    return list(visit([start]))


def oracle(graph, values, endpoints, thresholds):
    if any(values[v][k] < thresholds[k] for v in endpoints for k in range(len(thresholds))):
        return 'ineligible_endpoints'
    paths = simple_paths(graph,*endpoints)
    flags = [[all(values[v][k] >= thresholds[k] for v in path) for k in range(len(thresholds))] for path in paths]
    if any(all(row) for row in flags):
        return 'jointly_connected'
    if all(any(row[k] for row in flags) for k in range(len(thresholds))):
        return 'joint_only_obstruction'
    return 'single_trait_disconnection'


class SyntheticTests(unittest.TestCase):
    def test_all_16_patterns(self):
        positives = 0
        for a,b,c,d in itertools.product((0,1),repeat=4):
            values = [[1,1],[a,c],[b,d],[1,1]]
            result = rules.classify(SQUARE,values,(0,3),[1,1])['classification']
            direct = (a or b) and (c or d) and not ((a and c) or (b and d))
            self.assertEqual(result=='joint_only_obstruction',bool(direct))
            self.assertEqual(result,oracle(SQUARE,values,(0,3),[1,1]))
            positives += result=='joint_only_obstruction'
        self.assertEqual(positives,2)

    def test_all_ternary_two_trait_fields(self):
        #6561 complete fields, both endpoint diagonals, endpoint-relative bounds.
        for flat in itertools.product((-1,0,1),repeat=8):
            values = [list(flat[i:i+2]) for i in range(0,8,2)]
            for endpoints in ((0,3),(1,2)):
                tau = [min(values[endpoints[0]][k],values[endpoints[1]][k]) for k in range(2)]
                self.assertEqual(rules.classify(SQUARE,values,endpoints)['classification'],oracle(SQUARE,values,endpoints,tau))

    def test_additive_positive_and_aligned_negative(self):
        a=[1,0,2,1]; b=[1,2,0,1]
        self.assertEqual(rules.additive_projection(a)['contrast'],0)
        self.assertEqual(rules.additive_projection(b)['contrast'],0)
        self.assertEqual(rules.classify(SQUARE,list(map(list,zip(a,b))),(0,3))['classification'],'joint_only_obstruction')
        self.assertEqual(rules.classify(SQUARE,list(map(list,zip(a,a))),(0,3))['classification'],'jointly_connected')
        self.assertEqual(rules.classify(SQUARE,[[1,1],[0,1],[0,1],[1,1]],(0,3))['classification'],'single_trait_disconnection')

    def test_ties_and_fixed_ineligible_endpoints(self):
        self.assertEqual(rules.classify(SQUARE,[[1,1]]*4,(0,3))['classification'],'jointly_connected')
        self.assertEqual(rules.classify(SQUARE,[[1,1]]*4,(0,3),[2,1])['classification'],'ineligible_endpoints')
        values=[[1,1],[1-5e-11,2],[2,1-5e-11],[1,1]]
        self.assertEqual(rules.classify(SQUARE,values,(0,3))['classification'],'joint_only_obstruction')
        tol=rules.component_tolerances(values)
        self.assertEqual(tol,[2e-10,2e-10])
        self.assertEqual(rules.classify(SQUARE,values,(0,3),tolerance=tol)['classification'],'jointly_connected')

    def test_cube_bypass_and_unknown_vertices(self):
        values=[[1,1],[0,2],[2,0],[1,1]]+[[1,1]]*4
        self.assertEqual(rules.classify(SQUARE,values[:4],(0,3))['classification'],'joint_only_obstruction')
        self.assertEqual(rules.classify(CUBE,values,(0,3))['classification'],'jointly_connected')
        self.assertEqual(oracle(CUBE,values,(0,3),[1,1]),'jointly_connected')
        missing=values[:4]+[[None,None]]*4
        bound=rules.missing_vertex_bounds(CUBE,missing,(0,3))
        self.assertEqual(bound['status'],'inconclusive')
        self.assertTrue(bound['joint_optimistic']); self.assertFalse(bound['joint_pessimistic'])
        #An unknown trait cannot override a known failing trait at a vertex.
        partial=values[:4]+[[None,0]]*4
        self.assertEqual(rules.missing_vertex_bounds(CUBE,partial,(0,3))['status'],'certified_joint_only_obstruction')

    def test_interval_certificate_known_corners(self):
        robust=[[(1,2),(1,2)],[(0,.5),(1,2)],[(1,2),(0,.5)],[(1,2),(1,2)]]
        result=rules.interval_certificate(SQUARE,robust,(0,3))
        self.assertEqual(result['status'],'certified_joint_only_obstruction')
        uncertain=[row[:] for row in robust];uncertain[1]=[(0,1),(1,2)]
        self.assertEqual(rules.interval_certificate(SQUARE,uncertain,(0,3))['status'],'inconclusive')
        connected=[[(1,2),(1,2)]]*4
        self.assertEqual(rules.interval_certificate(SQUARE,connected,(0,3))['status'],'certified_jointly_connected')
        self.assertEqual(rules.interval_certificate(SQUARE,connected,(0,3),[1.5,1])['status'],'ineligible_definite_endpoints')
        for realization in itertools.product((0,1),repeat=8):
            field=[[robust[v][k][realization[2*v+k]] for k in range(2)] for v in range(4)]
            self.assertEqual(oracle(SQUARE,field,(0,3),result['thresholds']),'joint_only_obstruction')

    def test_interval_unknown_endpoint_and_invalid_input(self):
        with self.assertRaises(ValueError):rules.interval_certificate(SQUARE,[[None,(1,1)]]*4,(0,3))
        with self.assertRaises(ValueError):rules.interval_certificate(SQUARE,[[(2,1),(1,1)]]*4,(0,3))
        with self.assertRaises(ValueError):rules.classify(SQUARE,[[float('nan'),1]]*4,(0,3))
        with self.assertRaises(ValueError):rules.classify(SQUARE,[[1,1]]*4,(-1,3))
        with self.assertRaises(ValueError):rules.classify(SQUARE,[[1,1]]*4,(0,3),tolerance=-1)
        with self.assertRaises(ValueError):rules.classify(SQUARE,[[None,1]]*4,(0,3))

    def test_exact_rank_projection_and_saturation(self):
        design=np.array([[1,-1,-1],[1,-1,1],[1,1,-1],[1,1,1]],dtype=float)
        self.assertEqual(np.linalg.matrix_rank(design),3)
        self.assertEqual(4-np.linalg.matrix_rank(design),1)
        self.assertEqual(np.linalg.matrix_rank(design[:3]),3)
        self.assertEqual(3-np.linalg.matrix_rank(design[:3]),0)
        for values in ([Q(0),Q(0),Q(0),Q(4)],[Q(1,3),Q(-1,7),Q(5,2),Q(9,11)]):
            result=rules.additive_projection(values)
            self.assertEqual(result['projection'][0]-result['projection'][1]-result['projection'][2]+result['projection'][3],0)
            residual=np.array(list(map(float,result['residual'])))
            np.testing.assert_allclose(design.T@residual,0,atol=1e-14)
            coeff=np.linalg.lstsq(design,np.array(list(map(float,values))),rcond=None)[0]
            np.testing.assert_allclose(design@coeff,list(map(float,result['projection'])),atol=1e-14)
        self.assertEqual(rules.additive_projection([Q(0),Q(0),Q(0),Q(4)])['contrast'],4)

    def test_missing_corner_no_target_api(self):
        for intercept,a,b in itertools.product(range(-2,3),repeat=3):
            values=[Q(intercept),Q(intercept+b),Q(intercept+a),Q(intercept+a+b)]
            for missing in range(4):
                known={i:v for i,v in enumerate(values) if i!=missing}
                self.assertEqual(rules.missing_corner_predictions(known,missing)['additive'],values[missing])
        self.assertEqual(rules.missing_corner_predictions({0:10,1:1,2:1},3)['additive'],-8)
        with self.assertRaises(ValueError):rules.missing_corner_predictions({0:1,1:2,2:3,3:4},3)

    def test_threshold_cells_boundaries_and_unit_invariance(self):
        values=[[Q(1),Q(1)],[Q(0),Q(2)],[Q(2),Q(0)],[Q(1),Q(1)]]
        geometry=rules.threshold_geometry(SQUARE,values,(0,3))
        self.assertEqual(geometry['endpoint_caps'],[1,1])
        self.assertTrue(any(cell['classification']=='joint_only_obstruction' for cell in geometry['cells']))
        for cell in geometry['cells']:
            i,j=cell['axis_indices']
            tau=[geometry['axes'][0][i]['representative'],geometry['axes'][1][j]['representative']]
            self.assertEqual(cell['classification'],oracle(SQUARE,values,(0,3),tau))
        #Two strictly increasing transformations preserve endpoint-relative classes.
        transforms=[lambda v:v**3,lambda v:3*v+10]
        transformed=[[transforms[k](v) for k,v in enumerate(row)] for row in values]
        self.assertEqual(rules.classify(SQUARE,transformed,(0,3))['classification'],rules.classify(SQUARE,values,(0,3))['classification'])
        self.assertEqual(rules.native_shortfalls(values),{1:[1,0],2:[0,1]})

    def test_hash_convention_known_computation(self):
        ids=['a','b','c','d']
        integer=int(hashlib.sha256(b'C5-corner-v1|a|b|c|d').hexdigest(),16)
        self.assertEqual(rules.corner_index(ids),integer%4)
        with self.assertRaises(ValueError):rules.corner_index(['b','a','c','d'])


if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(SyntheticTests))
    receipt={'status':'passed' if result.wasSuccessful() else 'failed','tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
             'exact_reference_checks':{'binary_patterns':16,'ternary_field_diagonal_comparisons':13122,'interval_corner_realizations':256,'additive_missing_corner_recoveries':500},
             'data_scope':'Synthetic values only; no phenotype files opened','independence':'Separate simple-path and direct truth-table algorithms by same analyst, not blinded separate-person review',
             'hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('graph_rules.py','test_synthetic.py','metadata_split.py','split_manifest.json')}}
    (ROOT/'synthetic_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    raise SystemExit(0 if result.wasSuccessful() else 1)
