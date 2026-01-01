"""Checks only preselected sequence/design/source metadata, never assay outcomes."""
from pathlib import Path
import ast
import collections
import hashlib
import json


def main():
    folder = Path(__file__).resolve().parent
    csv_rows = json.loads((folder/'selected_library_metadata.json').read_text(encoding='utf-8'))
    xlsx_rows = json.loads((folder/'selected_workbook_metadata.json').read_text(encoding='utf-8'))
    cm = {x['library index']: x for x in csv_rows}
    if len(cm) != 9608 or len(xlsx_rows) != 9608:
        raise RuntimeError('Unexpected metadata cardinality')
    comparisons = {}
    for name in ['library index','subset','exonlength','exonstart in sequence','exonend in sequence','sequence']:
        comparisons[name] = sum(x.get(name) != cm[x['library index']].get(name) for x in xlsx_rows)
    if any(comparisons.values()):
        raise RuntimeError('Sequence/annotation mismatch between public source copies')
    gene_diffs = sum(x.get('gene name') != cm[x['library index']].get('gene name') for x in xlsx_rows)
    normalized_gene_diffs = sum((x.get('gene name') or '') != (cm[x['library index']].get('gene name') or '') for x in xlsx_rows)
    if normalized_gene_diffs:
        raise RuntimeError('Gene label mismatch beyond empty representation')
    source = ast.parse((folder/'author_code/ce.pipeline.py').read_text(encoding='utf-8'))
    constants = {}
    for node in ast.walk(source):
        if isinstance(node, ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
            name = node.targets[0].id
            if name in ['const0','const1','const2'] and isinstance(node.value,ast.Constant):
                constants[name] = node.value.value.upper()
    native = [r for r in xlsx_rows if r['subset']=='cassette_filtered']
    counts = {
        'rows':len(xlsx_rows), 'unique_library_IDs':len(cm),
        'sequence_lengths':dict(collections.Counter(len(x['sequence']) for x in xlsx_rows)),
        'unique_full210_sequences':len({x['sequence'].upper() for x in xlsx_rows}),
        'unique_barcodes_18_30':len({x['sequence'][18:30].upper() for x in xlsx_rows}),
        'unique_variable147_sequences_45_192':len({x['sequence'][45:192].upper() for x in xlsx_rows}),
        'nonempty_gene_labels':len({x['gene name'] for x in xlsx_rows if x.get('gene name')}),
        'unlabeled_gene_rows':sum(not x.get('gene name') for x in xlsx_rows),
        'subset_counts':dict(collections.Counter(x['subset'] for x in xlsx_rows)),
        'cassette_filtered_rows':len(native),
        'cassette_filtered_unique_genes':len({x['gene name'] for x in native}),
        'cassette_filtered_unique_variable_sequences':len({x['sequence'][45:192].upper() for x in native}),
        'cassette_filtered_gene_multiplicities':dict(collections.Counter(collections.Counter(x['gene name'] for x in native).values())),
        'exon_length_min':min(int(x['exonlength']) for x in xlsx_rows),
        'exon_length_max':max(int(x['exonlength']) for x in xlsx_rows),
        'exon_length_mod3_counts':dict(collections.Counter(int(x['exonlength'])%3 for x in xlsx_rows)),
        'exon_end_in_library_sequence_counts':dict(collections.Counter(x['exonend in sequence'] for x in xlsx_rows)),
    }
    if counts['exon_length_mod3_counts'] != {2:9608}:
        raise RuntimeError('Coordinate modulus check failed')
    map_info = {
        'coordinate_convention':'0-based half-open source slices; same convention checked in notebook source.',
        'source_to_reference_offset_after_library45':273,
        'amplicon_formula':'sequence[:30] + const1 + sequence[45:] + const2',
        'amplicon_length':30+len(constants['const1'])+165+len(constants['const2']),
        'const1_length':len(constants['const1']), 'const2_length':len(constants['const2']),
        'barcode_reference_interval':[18,30], 'variable_reference_interval':[318,465],
        'constant_upstream_donor':88, 'upstream_donor_dinucleotide':constants['const1'][58:60],
        'cassette_exon_start':'exonstart in sequence + 273', 'cassette_exon_end':435,
        'constant_downstream_acceptor':1214,
        'acceptor_preceding_dinucleotide':constants['const2'][729:731],
        'canonical_included_amplicon':'ref[:88] + ref[exonstart+273:435] + ref[1214:]',
        'canonical_skipped_amplicon':'ref[:88] + ref[1214:]',
        'full_ORF_status':'not_authenticated; amplication/mapping reference is not complete coding sequence',
        'unused_legacy_const0_length':len(constants['const0']),
        'do_not_infer':'Unused const0/tot_seq and legacy reporter tail are not proof of the 2026 promoter-specific full coding construct or terminating stop.',
    }
    for path in ['evidence.json','synthetic_selection_receipt.json']:
        if (folder/path).exists():
            raise RuntimeError('Refusing overwrite: '+path)
    obj = {'verdict':'NO-GO for authenticated full-CDS/synonymous accessibility; conditional metadata feasibility for matched-cassette RNA-context transfer only.',
           'sequence_metadata_counts':counts,'CSV_XLSX_exact_field_difference_counts':comparisons,
           'gene_name_empty_representation_differences':gene_diffs,'gene_name_differences_after_empty_normalization':normalized_gene_diffs,
           'intended_isoform_coordinate_map':map_info,
           'synonymous_edge_count':None,'synonymous_square_count':None,
           'graph_count_reason':'NOT COMPUTED: full coding phase, terminating stop and isoform protein consequence not authenticated. Null does not mean zero.',
           'assay_coverage_or_missingness':'Not accessed. Design identity is not measured context intersection.',
           'original_CV':'StratifiedKFold(n_splits=5,shuffle=True,random_state=42), no gene/context grouping passed at the inspected calls.',
           'protected_fields':'No assay phenotype, PSI, expression, outcome missingness or count data read.'}
    (folder/'evidence.json').write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    from selective_readers import self_test
    test = self_test()
    receipt = {'status':test,'selector_sha256':hashlib.sha256((folder/'selective_readers.py').read_bytes()).hexdigest(),
               'code_execution_scope':'Own selector and metadata checks only; deposited author code only AST-parsed/read, never executed.',
               'real_input_whitelist':['library index','subset','gene name','exonlength','exonstart in sequence','exonend in sequence','sequence']}
    (folder/'synthetic_selection_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'metadata_contract_checks_complete','counts':counts,'map':map_info,'evidence_sha256':hashlib.sha256((folder/'evidence.json').read_bytes()).hexdigest()}))


if __name__=='__main__':
    main()
