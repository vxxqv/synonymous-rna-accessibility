"""Invented-data checks only. Never opens the source workbook."""
import copy
import hashlib
import itertools
import json
import platform
import subprocess
import sys
import tempfile
import unittest
import zipfile
from fractions import Fraction
from pathlib import Path
from xml.sax.saxutils import escape

import analyze_three_states as a

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE.parent / "three_state_protocol.md"  # Text/protocol only.
N = a.NS["m"]
R = a.NS["r"]


def contract():
    return {
        "workbook_sha256": "", "sheet": "validation K562",
        "header_schema": {"C": "Paper_name",
            **{c: f"K562_gDNA_{i}_RPM" for i,c in enumerate(a.DNA, 1)},
            **{c: f"K562_cDNA_{i}_RPM" for i,c in enumerate(a.RNA, 1)}},
        "state_metadata": {
            key: {"row": a.ROWS[code], "label_cell": f"C{a.ROWS[code]}", "paper_label": a.LABELS[code],
                  "gDNA_RPM_cells": [f"{c}{a.ROWS[code]}" for c in a.DNA],
                  "cDNA_RPM_cells": [f"{c}{a.ROWS[code]}" for c in a.RNA]}
            for code,key in a.STATE_KEYS.items()
        },
    }


def numeric_inputs(amounts):
    records = {}
    for code in ("M", "E", "C"):
        dna, rna = amounts[code]
        for col,value in zip(a.DNA+a.RNA, list(dna)+list(rna)):
            cell = f"{col}{a.ROWS[code]}"
            records[cell] = a.numeric_record(cell, {"type":"n","lexical":str(value),
                                                    "has_formula":False,"inline":None})
    return records


def fixture(path, overrides=None, metadata_override=None, shared=False, sheet_name="validation K562"):
    cfg = contract()
    records = {}
    string_values = []
    def text_cell(address, text):
        if shared:
            index = len(string_values)
            string_values.append(text)
            return f'<c r="{address}" t="s"><v>{index}</v></c>'
        return f'<c r="{address}" t="inlineStr"><is><t>{escape(text)}</t></is></c>'
    for col,text in cfg["header_schema"].items():
        records[f"{col}1"] = text_cell(f"{col}1", text)
    for code,row in a.ROWS.items():
        records[f"C{row}"] = text_cell(f"C{row}", a.LABELS[code])
    for address in a.ALLOWED_NUMERIC:
        records[address] = f'<c r="{address}"><v>2</v></c>'
    records["L15"] = '<c r="L15" t="n"><v>DO_NOT_INSPECT_UNSELECTED_OUTCOME</v></c>'
    for address,xml in (overrides or {}).items():
        if xml is None:
            records.pop(address, None)
        else:
            records[address] = xml
    for address,text in (metadata_override or {}).items():
        records[address] = text_cell(address,text)
    rows = {}
    for address,xml in records.items():
        row = "".join(ch for ch in address if ch.isdigit())
        rows.setdefault(row,[]).append(xml)
    xml = f'<worksheet xmlns="{N}"><sheetData>'+''.join(
        f'<row r="{row}">{"".join(cells)}</row>' for row,cells in sorted(rows.items(),key=lambda x:int(x[0])))+'</sheetData></worksheet>'
    workbook = f'<workbook xmlns="{N}" xmlns:r="{R}"><sheets><sheet name="{sheet_name}" sheetId="1" r:id="rId1"/></sheets></workbook>'
    rels = '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml" Type="'+R+'/worksheet"/></Relationships>'
    with zipfile.ZipFile(path,"w") as z:
        z.writestr("xl/workbook.xml", workbook)
        z.writestr("xl/_rels/workbook.xml.rels",rels)
        z.writestr("xl/worksheets/sheet1.xml",xml)
        z.writestr("xl/worksheets/closed.xml","UNSELECTED_OUTCOME_SHEET_MUST_NOT_BE_PARSED")
        if shared:
            z.writestr("xl/sharedStrings.xml",f'<sst xmlns="{N}">'+''.join(
                f'<si><t>{escape(s)}</t></si>' for s in string_values)+'</sst>')
    cfg["workbook_sha256"] = a.sha(path)
    return cfg


def bounds(result, comparator):
    e = next(x["envelope"] for x in result["contrasts"] if x["contrast"]==comparator)
    return e


class Arithmetic(unittest.TestCase):
    def test_all_27_162_729_constant_fractions(self):
        inputs = numeric_inputs({"M":([2]*3,[6]*3),"E":([2]*3,[10]*3),"C":([2]*3,[2]*3)})
        result = a.evaluate(inputs)
        self.assertEqual(result["cardinalities"],{"states":3,"state_ratios":27,"contrast_quotients":162,"joint_shared_C_completions":729})
        for contrast,exact in [("C/M","1/3"),("C/E","1/5")]:
            b = bounds(result,contrast)
            self.assertEqual((b["lower"]["exact"],b["upper"]["exact"]),(exact,exact))
            self.assertTrue(b["strictly_lower_than_one"])
            self.assertTrue(b["at_least_twofold_lower"])
        self.assertEqual(result["joint_direction_counts"]["below_one|below_one"],729)
        self.assertEqual(sum(result["joint_direction_counts"].values()),729)

    def test_varied_independent_enumeration_every_case(self):
        amounts = {"M":([1,2,4],[2,4,8]), "E":([1,3,9],[3,9,27]), "C":([2,4,6],[1,2,3])}
        result = a.evaluate(numeric_inputs(amounts))
        reference = {
            code: {(i,j):Fraction(amounts[code][1][i-1],amounts[code][0][j-1])
                   for i,j in itertools.product(range(1,4),repeat=2)}
            for code in amounts
        }
        for row in result["state_results"]:
            for item in row["ratios"]:
                self.assertEqual(Fraction(item["value"]["exact"]),reference[row["state"]][item["RNA_index"],item["DNA_index"]])
        case_map = {r["id"]:r for r in result["contrast_quotients"]}
        for comp in ("M","E"):
            for i,j,k,l in itertools.product(range(1,4),repeat=4):
                key=f"C{comp}_{i}{j}_{k}{l}"
                expected=reference["C"][i,j]/reference[comp][k,l]
                self.assertEqual(Fraction(case_map[key]["value"]["exact"]),expected)
        for joint in result["joint_completions"]:
            self.assertEqual(case_map[joint["CM_case_id"]]["C_ratio_id"],joint["C_ratio_id"])
            self.assertEqual(case_map[joint["CE_case_id"]]["C_ratio_id"],joint["C_ratio_id"])
            self.assertEqual(case_map[joint["CM_case_id"]]["direction_at_one"],joint["CM_direction"])
            self.assertEqual(case_map[joint["CE_case_id"]]["direction_at_one"],joint["CE_direction"])
        self.assertEqual((bounds(result,"C/M")["lower"]["exact"],bounds(result,"C/M")["upper"]["exact"]),("1/48","3"))
        self.assertEqual((bounds(result,"C/E")["lower"]["exact"],bounds(result,"C/E")["upper"]["exact"]),("1/162","9/2"))
        self.assertTrue(all(d["envelope"]["strictly_lower_than_one"] for d in result["source_index_diagnostics"]))
        self.assertTrue(all(not d["agrees_with_primary_strict_flags"] for d in result["source_index_diagnostics"]))
        self.assertTrue(all(not r["envelope"]["strictly_lower_than_one"] for r in result["contrasts"]))

    def test_all_one_ties_are_not_strict(self):
        result=a.evaluate(numeric_inputs({s:([2]*3,[2]*3) for s in a.STATE_KEYS}))
        for row in result["contrasts"]:
            self.assertFalse(row["envelope"]["strictly_lower_than_one"])
            self.assertFalse(row["envelope"]["strictly_higher_than_one"])
            self.assertEqual(row["direction_case_counts"]["equal_one"],81)
            self.assertTrue(row["envelope"]["endpoint_equal_one"])

    def test_secondary_closed_boundaries(self):
        result=a.evaluate(numeric_inputs({"M":([1]*3,[2]*3),"E":([2]*3,[1]*3),"C":([1]*3,[1]*3)}))
        self.assertTrue(bounds(result,"C/M")["at_least_twofold_lower"])
        self.assertTrue(bounds(result,"C/E")["at_least_twofold_higher"])
        self.assertTrue(bounds(result,"C/M")["endpoint_equal_half_or_two"])
        self.assertTrue(bounds(result,"C/E")["endpoint_equal_half_or_two"])
        self.assertTrue(bounds(result,"C/M")["strictly_lower_than_one"])
        self.assertTrue(bounds(result,"C/E")["strictly_higher_than_one"])

    def test_zero_DNA_kept_unknown_no_dropout(self):
        result=a.evaluate(numeric_inputs({"M":([1]*3,[1]*3),"E":([1]*3,[1]*3),"C":([0,1,1],[1]*3)}))
        self.assertIn("UNKNOWN",result["status"])
        for row in result["contrasts"]:
            self.assertEqual(row["envelope"]["status"],"unknown")
            self.assertEqual(row["envelope"]["unknown_cases"],27)
            self.assertIsNone(row["envelope"]["strictly_lower_than_one"])
        self.assertEqual(len(result["joint_completions"]),729)

    def test_zero_RNA_is_recorded_zero(self):
        result=a.evaluate(numeric_inputs({"M":([1]*3,[1]*3),"E":([1]*3,[1]*3),"C":([1]*3,[0]*3)}))
        self.assertEqual(bounds(result,"C/M")["lower"]["exact"],"0")
        self.assertTrue(bounds(result,"C/M")["strictly_lower_than_one"])
        self.assertEqual(len(result["zero_source_cells"]),3)
        self.assertEqual(result["joint_direction_counts"]["below_one|below_one"],729)

    def test_positive_over_zero_and_zero_over_zero(self):
        positive=a.evaluate(numeric_inputs({"M":([1]*3,[0]*3),"E":([1]*3,[1]*3),"C":([1]*3,[1]*3)}))
        self.assertEqual(bounds(positive,"C/M")["lower"]["kind"],a.INF)
        self.assertTrue(bounds(positive,"C/M")["strictly_higher_than_one"])
        zero=a.evaluate(numeric_inputs({"M":([1]*3,[0]*3),"E":([1]*3,[1]*3),"C":([1]*3,[0]*3)}))
        self.assertEqual(bounds(zero,"C/M")["status"],"unknown")
        self.assertEqual(bounds(zero,"C/M")["unknown_cases"],81)
        self.assertTrue(bounds(zero,"C/E")["strictly_lower_than_one"])

    def test_one_bad_input_poison_primary_not_deleted(self):
        cells=numeric_inputs({s:([1]*3,[1]*3) for s in a.STATE_KEYS})
        cells["I15"]=a.numeric_record("I15",None)
        result=a.evaluate(cells)
        self.assertEqual(bounds(result,"C/M")["status"],"unknown")
        self.assertEqual(bounds(result,"C/M")["unknown_cases"],27)
        self.assertEqual(bounds(result,"C/E")["status"],"finite")
        self.assertEqual(len(result["source_cells"]),18)


class Selection(unittest.TestCase):
    def test_inline_and_shared_headers_exact_18_only(self):
        for shared in (False,True):
            with self.subTest(shared=shared),tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/"source.xlsx"
                cfg=fixture(path,shared=shared)
                seen=[]
                values,receipt=a.read_selected(path,cfg,seen.append)
                self.assertEqual(seen,list(a.ALLOWED_NUMERIC))
                self.assertEqual(set(values),set(a.ALLOWED_NUMERIC))
                self.assertFalse(receipt["unselected_numeric_cells_inspected"])

    def test_metadata_and_hash_guards_before_numeric(self):
        fixtures=[("header",{"F1":"WRONG"},"validation K562"),
                  ("label",{"C15":"WRONG"},"validation K562"),
                  ("sheet",{},"WRONG")]
        for name,meta,sheet in fixtures:
            with self.subTest(name=name),tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/"source.xlsx";cfg=fixture(path,metadata_override=meta,sheet_name=sheet);seen=[]
                with self.assertRaises(a.SchemaError):
                    a.read_selected(path,cfg,seen.append)
                self.assertEqual(seen,[])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"source.xlsx";cfg=fixture(path);cfg["workbook_sha256"]="0"*64;seen=[]
            with self.assertRaises(a.SchemaError):a.read_selected(path,cfg,seen.append)
            self.assertEqual(seen,[])

    def test_formula_numeric_cache_no_formula_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"source.xlsx"
            cfg=fixture(path,overrides={"F15":'<c r="F15"><f>CALL_UNTRUSTED_DO_NOT_EXECUTE()</f><v>2.50</v></c>'})
            values,_=a.read_selected(path,cfg)
            self.assertEqual(values["F15"]["exact"],"5/2")
            self.assertTrue(values["F15"]["has_formula"])
            self.assertNotIn("CALL_UNTRUSTED",json.dumps(values))

    def test_unknown_cache_missing_text_nonfinite_negative(self):
        examples={
            "uncached_formula":'<c r="F15"><f>2+3</f></c>',
            "absent_cell":None,
            "missing_numeric_value":'<c r="F15"/>',
            "nonnumeric_cell_type":'<c r="F15" t="str"><v>2</v></c>',
            "invalid_or_nonfinite_numeric_value":'<c r="F15"><v>NaN</v></c>',
            "negative_RPM":'<c r="F15"><v>-1</v></c>',
        }
        for reason,xml in examples.items():
            with self.subTest(reason=reason),tempfile.TemporaryDirectory() as tmp:
                path=Path(tmp)/"source.xlsx";cfg=fixture(path,overrides={"F15":xml})
                values,_=a.read_selected(path,cfg)
                self.assertEqual(values["F15"]["reason"],reason)
                result=a.evaluate(values)
                self.assertEqual(bounds(result,"C/M")["status"],"unknown")
                self.assertEqual(len(result["contrast_quotients"]),162)

    def test_decimal_lexical_exactness(self):
        for lexical,expected in [("0.1","1/10"),("1.23E-7","123/1000000000"),("9007199254740993","9007199254740993")]:
            row=a.numeric_record("F15",{"type":"n","lexical":lexical,"has_formula":False})
            self.assertEqual(row["exact"],expected)
            self.assertEqual(row["lexical"],lexical)
        invalid=a.numeric_record("F15",{"type":"n","lexical":"1/3","has_formula":False})
        self.assertEqual(invalid["kind"],a.UNKNOWN)

    def test_contract_whitelist_cannot_expand(self):
        cfg=contract()
        cfg["state_metadata"]["marker_only"]["gDNA_RPM_cells"].append("L15")
        with self.assertRaises(a.SchemaError):a.validate_contract(cfg)


class CLI(unittest.TestCase):
    def command(self,tmp,metadata_override=None):
        root=Path(tmp);workbook=root/"fixture.xlsx"
        cfg=fixture(workbook,metadata_override=metadata_override)
        protocol=root/"protocol.md";protocol.write_bytes(PROTOCOL.read_bytes())
        contract_path=root/"contract.json";a.dump_json(contract_path,cfg)
        freeze=root/"freeze.json"
        a.dump_json(freeze,{"implementation_sha256":a.sha(HERE/"analyze_three_states.py"),
                           "protocol_sha256":a.sha(protocol),"contract_sha256":a.sha(contract_path),
                           "workbook_sha256":cfg["workbook_sha256"],"runtime":{"python":platform.python_version()}})
        auth=root/"authorization.json"
        a.dump_json(auth,{"authorization_status":"APPROVED_EXACT_18_CELLS",
                         "freeze_sha256":a.sha(freeze),"authorized_numeric_cells":list(a.ALLOWED_NUMERIC),
                         "workbook_sha256":cfg["workbook_sha256"]})
        output=root/"out"
        cmd=[sys.executable,"-B",str(HERE/"analyze_three_states.py"),"--workbook",str(workbook),
             "--protocol",str(protocol),"--contract",str(contract_path),"--freeze",str(freeze),
             "--freeze-sha256",a.sha(freeze),"--authorization",str(auth),
             "--authorization-sha256",a.sha(auth),"--out",str(output)]
        return cmd,output,freeze

    def test_synthetic_cli_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            cmd,out,_=self.command(tmp)
            run=subprocess.run(cmd,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            self.assertEqual(run.stderr,"")
            result=a.load_json(out/"result.json")
            self.assertEqual(result["cardinalities"]["joint_shared_C_completions"],729)
            original={p.name:a.sha(p) for p in out.iterdir()}
            repeat=subprocess.run(cmd,capture_output=True,text=True)
            self.assertNotEqual(repeat.returncode,0)
            self.assertIn("no overwrite",repeat.stderr)
            self.assertEqual(original,{p.name:a.sha(p) for p in out.iterdir()})

    def test_optimized_refusal_before_output(self):
        for opt in ("-O","-OO"):
            with self.subTest(opt=opt),tempfile.TemporaryDirectory() as tmp:
                output=Path(tmp)/"must_not_exist"
                run=subprocess.run([sys.executable,opt,str(HERE/"analyze_three_states.py"),"--out",str(output)],
                                   capture_output=True,text=True)
                self.assertNotEqual(run.returncode,0)
                self.assertIn("Optimized execution",run.stderr)
                self.assertFalse(output.exists())

    def test_freeze_hash_refusal_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            cmd,out,_=self.command(tmp)
            cmd[cmd.index("--freeze-sha256")+1]="0"*64
            run=subprocess.run(cmd,capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0)
            self.assertIn("Freeze SHA256 mismatch",run.stderr)
            self.assertFalse(out.exists())

    def test_metadata_failure_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            cmd,out,_=self.command(tmp,{"C15":"WRONG"})
            run=subprocess.run(cmd,capture_output=True,text=True)
            self.assertNotEqual(run.returncode,0)
            failure=a.load_json(out/"failure.json")
            self.assertEqual(failure["status"],"FAILED_PRESERVED")
            self.assertFalse((out/"result.json").exists())
            self.assertIn("no numeric cells selected",failure["message"])


if __name__=="__main__":
    unittest.main(verbosity=2)
