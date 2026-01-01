"""Post-exposure synthetic parser/arithmetic checks, no source outcomes read."""
import sys
if not __debug__: raise RuntimeError('Optimization is refused')
sys.dont_write_bytecode=True
from decimal import Decimal,localcontext
import importlib.util
from io import StringIO
from pathlib import Path
import unittest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('source_audit',HERE/'audit_source.py')
a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)

class Tests(unittest.TestCase):
    def test_selected_cell_parser_does_not_decode_protected_cells(self):
        handler=a.Cells({2})
        a.sax_parse(StringIO('<worksheet><sheetData><row r="2"><c r="B2"><v>4</v></c><c r="E2"><v>PROTECTED</v></c><c r="L2"><v>PROTECTED</v></c></row><row r="3"><c r="B3"><v>UNAUTHORIZED</v></c></row></sheetData></worksheet>'),handler)
        self.assertEqual(set(handler.cells),{(2,'B')}); self.assertEqual(handler.cells[(2,'B')]['value_raw'],'4')
    def test_only_whitelisted_shared_strings_collected(self):
        handler=a.Strings({0})
        a.sax_parse(StringIO('<sst><si><t>ENSG1</t></si><si><t>PROTECTED</t></si></sst>'),handler)
        self.assertEqual(handler.values,{0:'ENSG1'})
    def test_formula_and_duplicate_selected_cells_refused(self):
        for xml in ('<x><c r="B2"><f>1+1</f><v>2</v></c></x>','<x><c r="B2"><v>2</v></c><c r="B2"><v>3</v></c></x>'):
            with self.assertRaises(ValueError): a.sax_parse(StringIO(xml),a.Cells({2}))
    def test_known_ratios_and_tenth_significant_half_even_ties(self):
        with localcontext() as context:
            context.prec=80
            self.assertEqual(a.logratio(Decimal(8),Decimal(4)),1)
            self.assertEqual(a.logratio(Decimal(3),Decimal(3)),0)
            self.assertEqual(a.significand10(Decimal('1.2345678905')),Decimal('1.234567890'))
            self.assertEqual(a.significand10(Decimal('1.2345678915')),Decimal('1.234567892'))
    def test_general_format_does_not_document_decimal_quantum(self):
        self.assertIsNone(a.display_quantum({'number_format':'General'}))
        self.assertEqual(a.display_quantum({'number_format':'0.000'}),Decimal('.001'))
        with localcontext() as context:
            context.prec=80
            self.assertTrue(a.interval_compatible(Decimal(10),Decimal(10),Decimal(0),Decimal('.01'),Decimal('.01'),Decimal('.001'))['compatible'])

if __name__=='__main__': unittest.main(verbosity=2)
