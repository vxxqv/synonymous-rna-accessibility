"""Synthetic tests only. No real workbook or outcome file is opened."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
import zipfile

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('reproduction',HERE/'reproduce.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def protocol():
    p=json.loads((HERE/'protocol.draft.json').read_text())
    p.update(expected_cohort_n=4,top_n=2,conventional_ks_pvalue=False)
    return p


def rows():
    return [dict(GENE_ID=gid,_row=i+2,HL_WT=1,HL_CNOT3KO=2**i,HL_FC_LOG2=i,
                 STAT_IP=3-i,LOG2_FC_IP=2*(3-i),SE_IP=2,
                 PTMD_SCORE=4-i,PTMD_SCORE_ABSOLUTE=2*i)
            for i,gid in enumerate(['A','B','C','D'])]


def workbook(path,data,poison=None,shared=False):
    ns=m.NS['s']; reg=m.NS['r']
    ET.register_namespace('',ns)
    root=ET.Element('{'+ns+'}worksheet'); sd=ET.SubElement(root,'{'+ns+'}sheetData')
    strings=[]
    def cell(parent,coord,value):
        if poison and coord in poison:
            c=ET.SubElement(parent,'{'+ns+'}c',r=coord)
            ET.SubElement(c,'{'+ns+'}f').text='POISON_UNREQUESTED_CELL()'
            ET.SubElement(c,'{'+ns+'}v').text='999'
        elif isinstance(value,str):
            if shared:
                strings.append(value)
                c=ET.SubElement(parent,'{'+ns+'}c',r=coord,t='s')
                ET.SubElement(c,'{'+ns+'}v').text=str(len(strings)-1)
            else:
                c=ET.SubElement(parent,'{'+ns+'}c',r=coord,t='inlineStr')
                ET.SubElement(ET.SubElement(c,'{'+ns+'}is'),'{'+ns+'}t').text=value
        elif value is not None:
            c=ET.SubElement(parent,'{'+ns+'}c',r=coord)
            ET.SubElement(c,'{'+ns+'}v').text=str(value)
    h=ET.SubElement(sd,'{'+ns+'}row',r='1')
    for col,name in m.FIELDS.items(): cell(h,col+'1',name)
    for i,record in enumerate(data,2):
        r=ET.SubElement(sd,'{'+ns+'}row',r=str(i))
        for col,name in m.FIELDS.items(): cell(r,col+str(i),record.get(name))
        cell(r,'L'+str(i),'UNREQUESTED_CONTEXT')
    wb=f'<workbook xmlns="{ns}" xmlns:r="{reg}"><sheets><sheet name="TableS2" sheetId="1" r:id="rId1"/></sheets></workbook>'
    rel='<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>'
    with zipfile.ZipFile(path,'w') as z:
        z.writestr('xl/workbook.xml',wb); z.writestr('xl/_rels/workbook.xml.rels',rel)
        z.writestr('xl/worksheets/sheet1.xml',ET.tostring(root))
        if shared:
            s=ET.Element('{'+ns+'}sst')
            for text in strings:
                ET.SubElement(ET.SubElement(s,'{'+ns+'}si'),'{'+ns+'}t').text=text
            z.writestr('xl/sharedStrings.xml',ET.tostring(s))


class Tests(unittest.TestCase):
    def test_known_cohorts_and_effects(self):
        a=m.analyze(rows(),protocol())
        self.assertEqual(a['high_ids'],['A','B']); self.assertEqual(a['n_cohort'],4)
        e=a['group_effects']['HL_FC_LOG2']
        self.assertEqual(e['median_difference'],-2); self.assertEqual(e['ks_statistic'],1)
        self.assertEqual(e['rank_biserial'],-1)
        self.assertEqual(a['spearman']['HL_FC_LOG2']['STAT_IP'],-1)
        self.assertEqual(a['spearman']['HL_FC_LOG2']['PTMD_SCORE_ABSOLUTE'],1)
        self.assertEqual(a['PTMD_overlap']['count'],2)
        self.assertFalse(a['discrepancies'])

    def test_ratio_inconsistency_retained(self):
        r=rows(); r[0]['HL_FC_LOG2']=0.01
        a=m.analyze(r,protocol())
        self.assertEqual(a['n_cohort'],4)
        self.assertIn('stored_HL_ratio_inconsistent:A',a['discrepancies'])

    def test_stat_not_silently_replaced(self):
        r=rows(); r[0]['LOG2_FC_IP']=-100
        a=m.analyze(r,protocol())
        self.assertEqual(a['high_ids'],['A','B'])
        self.assertIn('archived_IP_statistic_discrepancy:A',a['discrepancies'])

    def test_cutoff_ties(self):
        r=rows(); r[2]['STAT_IP']=2
        a=m.analyze(list(reversed(r)),protocol())
        self.assertEqual(a['high_ids'],['A','B'])
        self.assertTrue(a['selection_ties']['boundary_tie_crosses'])
        self.assertEqual(a['selection_ties']['all_cutoff_tie_ids'],['B','C'])

    def test_duplicates_not_deduplicated(self):
        r=rows(); r[1]['GENE_ID']='A'
        a=m.analyze(r,protocol())
        self.assertEqual(a['n_cohort'],2)
        self.assertIn('duplicate_GeneIDs_not_resolved',a['discrepancies'])
        self.assertEqual(sum('duplicate_gene_id' in x['reasons'] for x in a['exclusions']),2)

    def test_invalid_fields(self):
        for field,value in [('HL_WT',0),('HL_CNOT3KO',-1),('HL_WT',float('nan')),
                            ('HL_FC_LOG2',float('inf')),('STAT_IP',None),('LOG2_FC_IP','bad'),('SE_IP',float('inf'))]:
            with self.subTest(field=field,value=value):
                r=rows(); r[0][field]=value
                a=m.analyze(r,protocol())
                self.assertEqual(a['n_cohort'],3)
                self.assertTrue(any(x['GENE_ID']=='A' for x in a['exclusions']))

    def test_zero_negative_se_flagged_not_selected_away(self):
        for se in [0,-1]:
            r=rows(); r[0]['SE_IP']=se
            a=m.analyze(r,protocol())
            self.assertEqual(a['n_cohort'],4)
            self.assertTrue(any('SE_IP:A' in x for x in a['discrepancies']))

    def test_nonfinite_statistic_diagnostic_is_recorded(self):
        r=rows(); r[0]['LOG2_FC_IP']=1e308; r[0]['SE_IP']=1e-300
        a=m.analyze(r,protocol())
        self.assertIn('nonfinite_LFC_over_SE:A',a['discrepancies'])
        json.dumps(a,allow_nan=False)

    def test_missing_ptmd_no_alternative_overlap(self):
        r=rows(); r[0]['PTMD_SCORE']=None
        a=m.analyze(r,protocol())
        self.assertEqual(a['n_cohort'],4); self.assertEqual(a['spearman_complete_case_n'],3)
        self.assertEqual(a['PTMD_overlap']['status'],'not_computed')

    def test_rank_and_ks_known_answers(self):
        self.assertEqual(m.ranks([3,1,1,4]),[3,1.5,1.5,4])
        self.assertIsNone(m.spearman([1,1],[2,3]))
        self.assertEqual(m.ks_statistic([1,2],[2,3]),.5)
        self.assertEqual(m.distribution([1,2],[2,3])['probability_high_greater_with_half_ties'],.125)

    def test_log_ratio_handles_large_scale(self):
        r=rows()
        for x in r:
            x['HL_WT']*=1e200; x['HL_CNOT3KO']*=1e200
        self.assertFalse(any('HL_ratio' in x for x in m.analyze(r,protocol())['discrepancies']))

    def test_selected_xml_shared_and_inline(self):
        for shared in [False,True]:
            with self.subTest(shared=shared), tempfile.TemporaryDirectory() as td:
                p=Path(td)/'toy.xlsx'; workbook(p,rows(),poison={'L2','D5'},shared=shared)
                a=m.selected_cells(p,{'A','K'},{1,2})
                self.assertEqual(a[2]['GENE_ID'],'A'); self.assertEqual(a[2]['STAT_IP'],'3')
                self.assertNotIn(5,a)
                with self.assertRaisesRegex(ValueError,'Selected formula'):
                    m.selected_cells(p,{'D'},{5})

    def test_end_to_end_staged_exposure(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); wb=p/'toy.xlsx'; pp=p/'protocol.json'; fp=p/'freeze.json'
            r=rows()+[dict(rows()[0],GENE_ID='E',STAT_IP=None)]+[dict(rows()[1],GENE_ID='F',HL_WT=0)]
            workbook(wb,r,poison={'B6','C6','D6','F7','G7','L2'},shared=True)
            pr=protocol(); pr['source_workbook_sha256']=m.sha(wb)
            pp.write_text(json.dumps(pr))
            fp.write_text(json.dumps({'status':'root_authorized_outcome_read','protocol_sha256':m.sha(pp),
                                     'implementation_sha256':m.sha(HERE/'reproduce.py'),'workbook_sha256':m.sha(wb)}))
            a=m.execute(wb,pp,fp,p/'out')
            self.assertEqual(a['n_cohort'],4)
            exposure=json.loads((p/'out/exposure.json').read_text())
            decoded=sum([x['decoded_cells'] for x in exposure['stages']],[])
            for cell in ['6:B','6:C','6:D','7:F','7:G','2:L']:
                self.assertNotIn(cell,decoded)
            with self.assertRaises(FileExistsError): m.execute(wb,pp,fp,p/'out')

    def test_hash_and_authorization_gate_before_reader(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td); wb=p/'toy.xlsx'; pp=p/'protocol.json'; fp=p/'freeze.json'
            workbook(wb,rows()); pp.write_text(json.dumps(protocol()))
            for freeze in [{'status':'draft'}, {'status':'root_authorized_outcome_read','protocol_sha256':'bad'}]:
                fp.write_text(json.dumps(freeze))
                with patch.object(m,'selected_cells',side_effect=AssertionError('reader must not run')):
                    with self.assertRaises(ValueError): m.execute(wb,pp,fp,p/'out')

    def test_source_count_is_never_forced(self):
        p=protocol(); p['expected_cohort_n']=5481
        a=m.analyze(rows(),p)
        self.assertEqual(a['n_cohort'],4)
        self.assertIn('published_cohort_size_mismatch',a['discrepancies'])

    def test_optional_existing_scipy_reference(self):
        if importlib.util.find_spec('scipy') is None:
            self.skipTest('existing scipy unavailable; no installation')
        from scipy.stats import spearmanr,ks_2samp
        for a,b in [([1,2,3,4],[4,3,2,1]),([1,1,2,4],[1,2,2,3])]:
            self.assertAlmostEqual(m.spearman(a,b),float(spearmanr(a,b).statistic))
            self.assertAlmostEqual(m.ks_statistic(a,b),float(ks_2samp(a,b).statistic))
        self.assertIsNotNone(m.distribution([1,2],[2,3],True)['ks_pvalue'])


if __name__=='__main__':
    unittest.main(verbosity=2)
