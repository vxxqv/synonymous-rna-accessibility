import sys
if not __debug__: raise RuntimeError('Optimization is refused')
sys.dont_write_bytecode=True
import importlib.util
from pathlib import Path
import unittest
from Bio.Data import CodonTable

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('features',HERE/'build_features.py')
f=importlib.util.module_from_spec(spec); spec.loader.exec_module(f)

def fixture(strand='+'):
    r={'gene_id':'ENSG1.2','transcript_id':'ENST1.3','tags':[], 'strand':strand,
       'features':{'exon':[(100,116)]}}
    if strand=='+': r['features'].update({'CDS':[(102,110)],'start_codon':[(102,104)],'stop_codon':[(111,113)],'UTR':[(100,101),(111,116)]})
    else: r['features'].update({'CDS':[(106,114)],'start_codon':[(112,114)],'stop_codon':[(103,105)],'UTR':[(115,116),(100,105)]})
    return r,'ENST1.3|ENSG1.2|x|x|x|x|17|UTR5:1-2|CDS:3-14|UTR3:15-17|','CCATGGCTGCCTAAGGA'

class Tests(unittest.TestCase):
    def test_standard_code_all64_independently(self):
        table=CodonTable.unambiguous_dna_by_id[1]
        self.assertEqual(len(f.SENSE),61); self.assertEqual(len(f.AMINO),20)
        for c,aa in f.CODE.items(): self.assertEqual(aa,'*' if c in table.stop_codons else table.forward_table[c])
    def test_known_fractions_both_strands(self):
        for strand in ('+','-'):
            values,protein,meta=f.validate(*fixture(strand))
            self.assertEqual(protein,'MAA'); self.assertEqual(meta['CDS_spliced_0based_halfopen'],[2,14])
            self.assertEqual(values['CDS_length_nt'],12); self.assertEqual(values['CDS_GC'],.5)
            self.assertEqual(values['codon_ATG'],1/3); self.assertEqual(values['aa_A'],2/3)
            self.assertEqual(values['UTR5_GC'],1); self.assertEqual(values['UTR3_GC'],2/3)
            self.assertAlmostEqual(sum(values['codon_'+c] for c in f.SENSE),1)
            self.assertAlmostEqual(sum(values['aa_'+a] for a in f.AMINO),1)
    def test_spliced_mapping_and_split_codons(self):
        self.assertEqual(f.project([(100,104),(200,205)],'+',[(103,104),(200,202)]),[(3,8)])
        self.assertEqual(f.project([(100,104),(200,205)],'-',[(102,104),(200,201)]),[(4,9)])
    def test_ambiguous_and_internal_stops_and_start(self):
        for cds in ('ATGNNNTAA','ATGTGATAA','CTGGCTTAA','ATGGCTGCT','ATGGCTTA'):
            with self.assertRaises(ValueError): f.numeric(cds,'','')
        with self.assertRaises(ValueError): f.numeric('ATGTAA','N','')
    def test_zero_utrs_are_missing_GC_not_zero_GC(self):
        v,aa=f.numeric('ATGTAA','','')
        self.assertEqual(aa,'M'); self.assertEqual(v['UTR5_length_nt'],0)
        self.assertIsNone(v['UTR5_GC']); self.assertIsNone(v['UTR3_GC'])
    def test_mismatch_partial_and_sec(self):
        for mutation in ('header','partial','sec'):
            r,h,s=fixture()
            if mutation=='header': h=h.replace('CDS:3-14','CDS:4-14')
            if mutation=='partial': r['tags']=['cds_end_NF']
            if mutation=='sec': r['features']['Selenocysteine']=[(105,107)]
            with self.assertRaises(ValueError): f.validate(r,h,s)
    def test_stable_id_does_not_collapse_PAR_copy(self):
        self.assertEqual(f.stable('ENSG000001.2'),'ENSG000001')
        self.assertEqual(f.stable('ENSG000001.2_PAR_Y'),'ENSG000001.2_PAR_Y')
    def test_tie_is_literal_versioned_id_order(self):
        self.assertEqual(min([(-99,'ENST1.9'),(-99,'ENST1.10')]),(-99,'ENST1.10'))
    def test_exon_outside_and_overlap_fail(self):
        with self.assertRaises(ValueError): f.project([(1,5)],'+',[(4,6)])
        with self.assertRaises(ValueError): f.project([(1,5),(5,7)],'+',[(1,2)])

if __name__=='__main__': unittest.main(verbosity=2)
