"""Frozen selective K562 source audit, using standard-library exact arithmetic."""
from __future__ import annotations

if not __debug__:
    raise RuntimeError("Optimized execution (-O/-OO) is not permitted")

import argparse
import hashlib
import itertools
import json
import platform
import posixpath
import re
import sys
import traceback
import xml.etree.ElementTree as ET
import zipfile
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

PROTOCOL_SHA = "93feed845877ecbe31ad1297505c13cc9c01fbcd2966d38eff10fb3d4a30e366"
STATE_KEYS = {"M": "marker_only", "E": "edit_only", "C": "combined_A"}
ROWS = {"M": 15, "E": 16, "C": 17}
LABELS = {"M": "R1207R_ATT(G/A)", "E": "R1207R_(A/G)TTG", "C": "R1207R_(A/G)TT(G/A)"}
DNA = ("F", "G", "H")
RNA = ("I", "J", "K")
ALLOWED_NUMERIC = tuple(f"{col}{row}" for row in (15, 16, 17) for col in DNA + RNA)
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
TAG = "{" + NS["m"] + "}"
FINITE, INF, UNKNOWN = "finite", "positive_infinity", "unknown"


class SchemaError(RuntimeError):
    pass


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def dump_json(path, data):
    path = Path(path)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


def finite(value):
    return {"kind": FINITE, "value": Fraction(value), "reason": None}


def unknown(reason):
    return {"kind": UNKNOWN, "value": None, "reason": reason}


def infinity(reason="positive_source_ratio_over_zero"):
    return {"kind": INF, "value": None, "reason": reason}


def decimal_text(value):
    with localcontext() as context:
        context.prec = 50
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def encoded(value):
    return {
        "kind": value["kind"],
        "exact": str(value["value"]) if value["kind"] == FINITE else "+inf" if value["kind"] == INF else None,
        "decimal_display_50digits": decimal_text(value["value"]) if value["kind"] == FINITE else None,
        "reason": value["reason"],
    }


def divide(top, bottom):
    if top["kind"] == UNKNOWN or bottom["kind"] == UNKNOWN:
        return unknown("unknown_state_ratio_in_contrast")
    if top["kind"] != FINITE or bottom["kind"] != FINITE:
        raise RuntimeError("State ratios must be finite or unknown")
    if bottom["value"] == 0:
        return infinity() if top["value"] > 0 else unknown("zero_over_zero")
    return finite(top["value"] / bottom["value"])


def direction(value):
    if value["kind"] == UNKNOWN:
        return "unknown"
    if value["kind"] == INF or value["value"] > 1:
        return "above_one"
    if value["value"] < 1:
        return "below_one"
    return "equal_one"


def envelope(values):
    if not values:
        raise RuntimeError("Empty completion set")
    missing = sum(v["kind"] == UNKNOWN for v in values)
    if missing:
        return {"status": "unknown", "lower": encoded(finite(0)), "upper": encoded(infinity("unbounded_unknown")),
                "unknown_cases": missing, "total_cases": len(values),
                "strictly_lower_than_one": None, "strictly_higher_than_one": None,
                "at_least_twofold_lower": None, "at_least_twofold_higher": None,
                "endpoint_equal_one": None, "endpoint_equal_half_or_two": None}
    present = [v["value"] for v in values if v["kind"] == FINITE]
    low = finite(min(present)) if present else infinity()
    high = infinity() if any(v["kind"] == INF for v in values) else finite(max(present))
    lower_inf, upper_inf = low["kind"] == INF, high["kind"] == INF
    return {
        "status": "defined_extended" if upper_inf else "finite",
        "lower": encoded(low), "upper": encoded(high), "unknown_cases": 0, "total_cases": len(values),
        "strictly_lower_than_one": not upper_inf and high["value"] < 1,
        "strictly_higher_than_one": lower_inf or low["value"] > 1,
        "at_least_twofold_lower": not upper_inf and high["value"] <= Fraction(1, 2),
        "at_least_twofold_higher": lower_inf or low["value"] >= 2,
        "endpoint_equal_one": any(v["kind"] == FINITE and v["value"] == 1 for v in (low, high)),
        "endpoint_equal_half_or_two": any(v["kind"] == FINITE and v["value"] in (Fraction(1, 2), Fraction(2)) for v in (low, high)),
    }


def validate_contract(contract):
    if contract["sheet"] != "validation K562" or set(contract["state_metadata"]) != set(STATE_KEYS.values()):
        raise SchemaError("Wrong worksheet or state inventory")
    required = {"C": "Paper_name", **{col: f"K562_gDNA_{i}_RPM" for i, col in enumerate(DNA, 1)},
                **{col: f"K562_cDNA_{i}_RPM" for i, col in enumerate(RNA, 1)}}
    if contract["header_schema"] != required:
        raise SchemaError("Wrong frozen header contract")
    selected = []
    for code, key in STATE_KEYS.items():
        state = contract["state_metadata"][key]
        row = ROWS[code]
        if (state["row"], state["label_cell"], state["paper_label"]) != (row, f"C{row}", LABELS[code]):
            raise SchemaError("Wrong state/row/label contract")
        if state["gDNA_RPM_cells"] != [f"{c}{row}" for c in DNA] or state["cDNA_RPM_cells"] != [f"{c}{row}" for c in RNA]:
            raise SchemaError("Numeric whitelist differs from approved cells")
        selected.extend(state["gDNA_RPM_cells"] + state["cDNA_RPM_cells"])
    if set(selected) != set(ALLOWED_NUMERIC) or len(selected) != 18:
        raise SchemaError("Not exactly the 18 approved numeric cells")


def worksheet_member(archive, sheet_name):
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    sheets = [s for s in workbook.findall("m:sheets/m:sheet", NS) if s.get("name") == sheet_name]
    if len(sheets) != 1:
        raise SchemaError("Worksheet name missing or duplicate")
    rid = sheets[0].get("{" + NS["r"] + "}id")
    rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    matches = [r for r in rels if r.get("Id") == rid]
    if len(matches) != 1 or matches[0].get("TargetMode") == "External":
        raise SchemaError("Worksheet relationship missing, duplicate or external")
    target = matches[0].get("Target", "")
    member = posixpath.normpath(target.lstrip("/") if target.startswith("/") else posixpath.join("xl", target))
    if not member.startswith("xl/") or member not in archive.namelist():
        raise SchemaError("Invalid worksheet member")
    return member


def select_cells(archive, member, wanted):
    """Access v/is/f children only for the selected addresses."""
    found = {}
    with archive.open(member) as stream:
        for event, cell in ET.iterparse(stream, events=("end",)):
            if cell.tag != TAG + "c":
                continue
            address = cell.get("r")
            if address in wanted:
                if address in found:
                    raise SchemaError(f"Duplicate selected cell {address}")
                v = cell.find(TAG + "v")
                inline = cell.find(TAG + "is")
                found[address] = {
                    "type": cell.get("t", "n"),
                    "lexical": None if v is None else v.text,
                    "inline": None if inline is None else "".join(inline.itertext()),
                    "has_formula": cell.find(TAG + "f") is not None,
                }
            cell.clear()
    return found


def metadata_strings(archive, cells):
    wanted = set()
    for cell in cells.values():
        if cell["has_formula"]:
            raise SchemaError("Formula in an identity/header cell")
        if cell["type"] == "s":
            try:
                wanted.add(int(cell["lexical"]))
            except (TypeError, ValueError):
                raise SchemaError("Invalid shared-string metadata index") from None
    strings = {}
    if wanted:
        if "xl/sharedStrings.xml" not in archive.namelist():
            raise SchemaError("Missing shared strings")
        with archive.open("xl/sharedStrings.xml") as stream:
            index = 0
            for event, node in ET.iterparse(stream, events=("end",)):
                if node.tag == TAG + "si":
                    if index in wanted:
                        strings[index] = "".join(t.text or "" for t in node.iter(TAG + "t"))
                    index += 1
                    node.clear()
        if set(strings) != wanted:
            raise SchemaError("Missing shared-string metadata entry")
    result = {}
    for address, cell in cells.items():
        if cell["type"] == "s":
            result[address] = strings[int(cell["lexical"])]
        elif cell["type"] == "inlineStr":
            result[address] = cell["inline"]
        elif cell["type"] == "str":
            result[address] = cell["lexical"]
        else:
            raise SchemaError("Nonstring identity/header cell")
    return result


def numeric_record(address, cell):
    if cell is None:
        return {"cell": address, "kind": UNKNOWN, "lexical": None, "exact": None,
                "has_formula": False, "reason": "absent_cell"}
    record = {"cell": address, "type": cell["type"], "lexical": cell["lexical"],
              "has_formula": cell["has_formula"], "kind": UNKNOWN, "exact": None, "reason": None}
    if cell["has_formula"] and cell["lexical"] is None:
        record["reason"] = "uncached_formula"
    elif cell["type"] != "n":
        record["reason"] = "nonnumeric_cell_type"
    elif cell["lexical"] is None:
        record["reason"] = "missing_numeric_value"
    elif re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?", cell["lexical"].strip()) is None:
        record["reason"] = "invalid_or_nonfinite_numeric_value"
    else:
        try:
            amount = Fraction(cell["lexical"])
        except (ValueError, ZeroDivisionError, TypeError):
            record["reason"] = "invalid_or_nonfinite_numeric_value"
        else:
            if amount < 0:
                record["reason"] = "negative_RPM"
            else:
                record.update(kind=FINITE, exact=str(amount))
    return record


def read_selected(workbook_path, contract, on_numeric_read=None):
    validate_contract(contract)
    if sha(workbook_path) != contract["workbook_sha256"]:
        raise SchemaError("Workbook SHA256 mismatch")
    with zipfile.ZipFile(workbook_path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise SchemaError("Duplicate ZIP members")
        member = worksheet_member(archive, contract["sheet"])
        expected = {f"{col}1": name for col, name in contract["header_schema"].items()}
        expected.update({f"C{ROWS[code]}": LABELS[code] for code in STATE_KEYS})
        raw_metadata = select_cells(archive, member, set(expected))
        if set(raw_metadata) != set(expected):
            raise SchemaError("Missing metadata cells")
        if metadata_strings(archive, raw_metadata) != expected:
            raise SchemaError("Source labels/header mismatch; no numeric cells selected")
        # This second pass is entered only after every identity/header guard.
        raw_numbers = select_cells(archive, member, set(ALLOWED_NUMERIC))
        values = {}
        for address in ALLOWED_NUMERIC:
            if on_numeric_read:
                on_numeric_read(address)
            values[address] = numeric_record(address, raw_numbers.get(address))
    return values, {"worksheet_member": member, "metadata_cells_verified": sorted(expected),
                    "numeric_cells_selected": list(ALLOWED_NUMERIC), "numeric_count": 18,
                    "unselected_numeric_cells_inspected": False}


def source_value(record):
    return finite(Fraction(record["exact"])) if record["kind"] == FINITE else unknown(record["reason"])


def evaluate(cells):
    if set(cells) != set(ALLOWED_NUMERIC):
        raise SchemaError("Numeric record keys differ from fixed whitelist")
    ratios = {}
    state_rows = []
    for code in ("M", "E", "C"):
        row = ROWS[code]
        entries = []
        for i, j in itertools.product(range(1, 4), repeat=2):
            rcell, dcell = f"{RNA[i-1]}{row}", f"{DNA[j-1]}{row}"
            rv, dv = source_value(cells[rcell]), source_value(cells[dcell])
            if rv["kind"] == UNKNOWN or dv["kind"] == UNKNOWN:
                value = unknown("unknown_source_RPM")
            elif dv["value"] == 0:
                value = unknown("zero_DNA_RPM")
            else:
                value = finite(rv["value"] / dv["value"])
            key = f"{code}_R{i}_D{j}"
            ratios[key] = value
            entries.append({"id": key, "state": code, "RNA_index": i, "DNA_index": j,
                            "RNA_cell": rcell, "DNA_cell": dcell, "value": encoded(value)})
        state_rows.append({"state": code, "ratios": entries,
                           "envelope": envelope([ratios[e["id"]] for e in entries]),
                           "source_index_ids": [f"{code}_R{i}_D{i}" for i in range(1, 4)]})

    cases = []
    lookup = {}
    contrast_rows = []
    diagnostics = []
    for comparator in ("M", "E"):
        ids, diagonal_ids = [], []
        values, diagonal_values = [], []
        for i, j, k, l in itertools.product(range(1, 4), repeat=4):
            top, bottom = f"C_R{i}_D{j}", f"{comparator}_R{k}_D{l}"
            key = f"C{comparator}_{i}{j}_{k}{l}"
            value = divide(ratios[top], ratios[bottom])
            lookup[(comparator, top, bottom)] = (key, value)
            entry = {"id": key, "contrast": f"C/{comparator}", "C_ratio_id": top,
                     "comparator_ratio_id": bottom, "value": encoded(value),
                     "direction_at_one": direction(value)}
            cases.append(entry)
            ids.append(key)
            values.append(value)
            if i == j and k == l:
                diagonal_ids.append(key)
                diagonal_values.append(value)
        primary_envelope = envelope(values)
        counts = {name: sum(direction(v) == name for v in values)
                  for name in ("below_one", "equal_one", "above_one", "unknown")}
        contrast_rows.append({"contrast": f"C/{comparator}", "case_ids": ids,
                              "envelope": primary_envelope, "direction_case_counts": counts})
        diagonal_envelope = envelope(diagonal_values)
        diagnostics.append({
            "contrast": f"C/{comparator}", "scope": "index-matched within-state ratios; no cross-state pairing assumed",
            "case_ids": diagonal_ids, "envelope": diagonal_envelope,
            "agrees_with_primary_strict_flags": (
                diagonal_envelope["strictly_lower_than_one"] == primary_envelope["strictly_lower_than_one"]
                and diagonal_envelope["strictly_higher_than_one"] == primary_envelope["strictly_higher_than_one"]),
        })

    joint = []
    directions = ("below_one", "equal_one", "above_one", "unknown")
    joint_counts = {f"{a}|{b}": 0 for a, b in itertools.product(directions, repeat=2)}
    for ci, cj, mi, mj, ei, ej in itertools.product(range(1, 4), repeat=6):
        ckey, mkey, ekey = f"C_R{ci}_D{cj}", f"M_R{mi}_D{mj}", f"E_R{ei}_D{ej}"
        cm_id, cm = lookup[("M", ckey, mkey)]
        ce_id, ce = lookup[("E", ckey, ekey)]
        dcm, dce = direction(cm), direction(ce)
        pair = f"{dcm}|{dce}"
        joint_counts[pair] = joint_counts.get(pair, 0) + 1
        joint.append({"id": f"J_{ci}{cj}_{mi}{mj}_{ei}{ej}", "C_ratio_id": ckey, "M_ratio_id": mkey,
                      "E_ratio_id": ekey, "CM_case_id": cm_id, "CE_case_id": ce_id,
                      "CM_direction": dcm, "CE_direction": dce})
    cardinalities = {"states": len(state_rows), "state_ratios": len(ratios),
                     "contrast_quotients": len(cases), "joint_shared_C_completions": len(joint)}
    if cardinalities != {"states": 3, "state_ratios": 27, "contrast_quotients": 162, "joint_shared_C_completions": 729}:
        raise RuntimeError("Completion cardinality mismatch")
    unknown_inputs = [r for r in cells.values() if r["kind"] == UNKNOWN]
    zero_inputs = [r["cell"] for r in cells.values() if r["kind"] == FINITE and Fraction(r["exact"]) == 0]
    has_unknown = bool(unknown_inputs) or any(v["kind"] == UNKNOWN for v in ratios.values()) or any(
        r["envelope"]["status"] == "unknown" for r in contrast_rows)
    return {
        "status": "COMPLETED_WITH_UNKNOWN_MEASUREMENTS_OR_RATIOS" if has_unknown else "COMPLETED",
        "cardinalities": cardinalities,
        "source_cells": [cells[c] for c in ALLOWED_NUMERIC],
        "unknown_source_cells": unknown_inputs, "zero_source_cells": zero_inputs,
        "state_results": state_rows, "contrasts": contrast_rows, "contrast_quotients": cases,
        "source_index_diagnostics": diagnostics, "joint_completions": joint,
        "joint_direction_counts": joint_counts,
        "scope": [
            "Finite-source conditional-background ratios, not isolated intended-edit effects or interactions",
            "Full Cartesian pairing does not establish culture independence or remove background/recovery differences",
            "Empirical envelopes are not confidence intervals; completions are not independent replicates",
            "Zero/infinity arithmetic describes recorded assay values, not biological zero/infinite effects",
            "Original qualitative example was published and exposed before this numerical audit",
        ],
    }


def verify_freeze(args):
    if sha(args.freeze) != args.freeze_sha256:
        raise SchemaError("Freeze SHA256 mismatch")
    freeze = load_json(args.freeze)
    if sha(__file__) != freeze["implementation_sha256"]:
        raise SchemaError("Implementation SHA256 mismatch")
    if sha(args.protocol) != freeze["protocol_sha256"] or freeze["protocol_sha256"] != PROTOCOL_SHA:
        raise SchemaError("Protocol SHA256 mismatch")
    if sha(args.contract) != freeze["contract_sha256"]:
        raise SchemaError("Contract SHA256 mismatch")
    if platform.python_implementation() != "CPython" or platform.python_version() != freeze["runtime"]["python"]:
        raise SchemaError("Runtime mismatch")
    authorization = load_json(args.authorization)
    if sha(args.authorization) != args.authorization_sha256:
        raise SchemaError("Authorization receipt SHA256 mismatch")
    if authorization["freeze_sha256"] != args.freeze_sha256 or authorization["authorized_numeric_cells"] != list(ALLOWED_NUMERIC):
        raise SchemaError("Authorization scope mismatch")
    if authorization["authorization_status"] != "APPROVED_EXACT_18_CELLS":
        raise SchemaError("Authorization is not approved")
    if authorization["workbook_sha256"] != freeze["workbook_sha256"]:
        raise SchemaError("Authorization source mismatch")
    contract = load_json(args.contract)
    if contract["workbook_sha256"] != freeze["workbook_sha256"]:
        raise SchemaError("Frozen contract source mismatch")
    validate_contract(contract)
    return freeze, contract


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workbook", "protocol", "contract", "freeze", "freeze-sha256",
                 "authorization", "authorization-sha256", "out"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    output = Path(args.out)
    if output.exists():
        raise FileExistsError("Output path already exists; no overwrite permitted")
    freeze, contract = verify_freeze(args)
    output.mkdir(parents=False, exist_ok=False)
    try:
        values, selection = read_selected(args.workbook, contract)
        result = evaluate(values)
        if sha(args.workbook) != contract["workbook_sha256"]:
            raise SchemaError("Source workbook changed after selection")
        result.update(protocol_sha256=freeze["protocol_sha256"], contract_sha256=freeze["contract_sha256"],
                      freeze_sha256=args.freeze_sha256, source_workbook_sha256=sha(args.workbook),
                      implementation_sha256=sha(__file__), selection=selection,
                      runtime={"implementation": platform.python_implementation(), "python": platform.python_version(),
                               "platform": platform.platform()})
        dump_json(output / "result.json", result)
        dump_json(output / "receipt.json", {
            "status": result["status"], "result_sha256": sha(output / "result.json"),
            "source_sha256_before_after": sha(args.workbook), "freeze_sha256": args.freeze_sha256,
            "authorization_sha256": args.authorization_sha256,
            "selected_numeric_cells": list(ALLOWED_NUMERIC), "cardinalities": result["cardinalities"],
            "new_normalization_or_pseudocount": False, "source_score_cells_opened": False,
        })
        dump_json(output / "manifest.json", {
            "files": {n: sha(output / n) for n in ("result.json", "receipt.json")},
            "status": result["status"],
        })
        print(json.dumps({"status": result["status"], "result_sha256": sha(output / "result.json"),
                          "cardinalities": result["cardinalities"]}, indent=2))
    except Exception as exc:
        dump_json(output / "failure.json", {"status": "FAILED_PRESERVED", "type": type(exc).__name__,
                                           "message": str(exc), "traceback": traceback.format_exc()})
        raise


if __name__ == "__main__":
    main()
