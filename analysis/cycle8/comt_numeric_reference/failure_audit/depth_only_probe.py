"""Frozen failure-only probe: decode identity/DP and RF status fields exclusively."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited before input access')

import argparse
from collections import defaultdict, Counter
import csv
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path

EXTRACT_SEAL = '4e365b305c356a25762c540092083c1a6e87bbc89041c97eb90a2e8f0436d7f0'
RF_SEAL = '091b66e102669f466121244e65a85c5434eafd514e903285074515a51bdd7564'
CONTRACT_SHA = '091ab6a57b7d64ab56638430e95befca5dbffd790eec9c62ea5d6a504d450ce5'
IDENTITY_DP = ('record_id', 'Sample', '#CHROM', 'POS', 'REF', 'ALT', 'POS_NT', 'REF_NT', 'ALT_NT', 'DP')
RF_FIELDS = ('record_id', 'VAR_ID', 'frequency_pass', 'globally_rejected')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def selected_rows(stream, names):
    # Binary split only. Unselected numeric-field bytes are never decoded or parsed.
    header = stream.readline().rstrip(b'\r\n').split(b'\t')
    if len(set(header)) != len(header):
        raise ValueError('Duplicate header')
    positions = {}
    for name in names:
        key = name.encode('ascii')
        if key not in header:
            raise ValueError('Required selected field missing: ' + name)
        positions[name] = header.index(key)
    for line in stream:
        fields = line.rstrip(b'\r\n').split(b'\t')
        if len(fields) != len(header):
            raise ValueError('Malformed TSV structure')
        yield {name: fields[i].decode('utf8') for name, i in positions.items()}


def med(x):
    s = sorted(x)
    return s[len(s) // 2] if len(s) % 2 else (s[len(s) // 2 - 1] + s[len(s) // 2]) / 2


def key_for_target(first, perbase):
    if ',' in first['REF'] or ',' in first['ALT']:
        pairs = sorted({(int(r['POS_NT']), r['REF_NT'], r['ALT_NT']) for r in perbase if r['REF_NT'] != r['ALT_NT']})
    else:
        if len(first['REF']) != len(first['ALT']):
            raise ValueError('Unsupported unequal-length source representation')
        pairs = [(int(first['POS']) + i, a, b) for i, (a, b) in enumerate(zip(first['REF'], first['ALT'])) if a != b]
    return ';'.join(str(p) + ':' + a + '>' + b for p, a, b in pairs)


def reconstruct(row_iterables, mask_iterables, contract):
    masks = {}
    for it in mask_iterables:
        for r in it:
            if r['record_id'] in masks:
                raise ValueError('Duplicate RF record ID')
            if r['frequency_pass'] not in ('TRUE', 'FALSE') or r['globally_rejected'] not in ('TRUE', 'FALSE'):
                raise ValueError('Unexpected mask flag')
            masks[r['record_id']] = r
    members = {v['sample']: v for v in contract['members'].values()}
    if len(members) != len(contract['members']):
        raise ValueError('Duplicate sample metadata')
    targets = {(int(r['ROI']), r['edit_key']): r for r in contract['targets']}
    if len(targets) != len(contract['targets']):
        raise ValueError('Duplicate target key')
    by_allele, seen = defaultdict(list), set()
    counts = Counter()
    for it in row_iterables:
        for r in it:
            counts['selected_extracted_rows'] += 1
            rid = r['record_id']
            if rid in seen or rid not in masks:
                raise ValueError('Duplicate row ID or unmasked row')
            seen.add(rid)
            vid = ':'.join(r[k] for k in ('#CHROM', 'POS', 'REF', 'ALT'))
            if masks[rid]['VAR_ID'] != vid:
                raise ValueError('Mask/source raw identity mismatch')
            if r['Sample'] not in members:
                raise ValueError('Sample missing from sealed metadata')
            if masks[rid]['frequency_pass'] == 'TRUE' and masks[rid]['globally_rejected'] == 'FALSE':
                r['DP'] = F(r['DP'])
                if r['DP'] < 0:
                    raise ValueError('Negative depth')
                by_allele[(r['Sample'], vid)].append(r)
                counts['retained_per_base_rows'] += 1
    if seen != set(masks):
        raise ValueError('RF masks and extracted record universe differ')
    groups = defaultdict(list)
    for (sample, vid), rows in by_allele.items():
        first = rows[0]
        depth = med([r['DP'] for r in rows])
        roi = int(members[sample]['ROI'])
        target = targets.get((roi, key_for_target(first, rows)))
        # SNP here means one-base raw REF/ALT. Multi-base one-change records remain MNPs.
        raw_snv = len(first['REF']) == len(first['ALT']) == 1 and ',' not in first['ALT']
        groups[(sample, int(first['POS']))].append({'VAR_ID': vid, 'DP': depth,
            'per_base_depth_values': sorted({r['DP'] for r in rows}),
            'raw_SNV': raw_snv, 'target': target['hgvs_nt'] if target else None,
            'target_SNV': bool(target['is_SNV']) if target else None,
            'trait': members[sample]['trait'], 'ROI': roi})
        counts['collapsed_retained_alleles'] += 1
        counts['collapsed_alleles_with_per_base_DP_variation'] += len({r['DP'] for r in rows}) > 1
    conflicts = {k: v for k, v in groups.items() if len({a['DP'] for a in v}) > 1}
    impacts_all, impacts_assay = {}, {}
    snv_consistent = lower_mnp = first_snv = 0
    table = []
    for (sample, pos), alleles in sorted(conflicts.items()):
        dps = sorted({a['DP'] for a in alleles})
        snv = {a['DP'] for a in alleles if a['raw_SNV']}
        mnp = {a['DP'] for a in alleles if not a['raw_SNV']}
        snv_consistent += len(snv) == 1
        lower_mnp += bool(snv) and bool(mnp) and max(mnp) < min(snv)
        first = min(alleles, key=lambda a: a['VAR_ID'])
        first_snv += first['raw_SNV']
        target_names = sorted({a['target'] for a in alleles if a['target'] is not None})
        for a in alleles:
            if a['target'] is not None:
                rec = {'hgvs_nt': a['target'], 'Sample': sample, 'POS': pos, 'trait': a['trait'],
                       'target_SNV': a['target_SNV'], 'raw_id': a['VAR_ID'], 'allele_DP': str(a['DP'])}
                impacts_all[(a['target'], sample)] = rec
                if a['trait'] != 'NC':
                    impacts_assay[(a['target'], sample)] = rec
        table.append({'Sample': sample, 'POS': pos, 'trait': alleles[0]['trait'],
            'unique_collapsed_depths': ';'.join(str(d) for d in dps),
            'collapsed_alleles': len(alleles), 'target_count': len(target_names),
            'target_ids': ';'.join(target_names), 'source_lexical_first_VAR_ID': first['VAR_ID'],
            'source_lexical_first_DP': str(first['DP']), 'source_lexical_first_is_raw_SNV': first['raw_SNV'],
            'all_depth_and_identity_records': json.dumps(sorted(alleles, key=lambda a: a['VAR_ID']), default=str)})
    def coverage(records):
        return {'target_sample_records': len(records), 'target_IDs': len({v['hgvs_nt'] for v in records.values()}),
            'SNV_target_IDs': len({v['hgvs_nt'] for v in records.values() if v['target_SNV']}),
            'nonSNV_target_IDs': len({v['hgvs_nt'] for v in records.values() if not v['target_SNV']})}
    summary = {'counts': dict(counts), 'retained_sample_start_groups': len(groups),
        'conflicting_groups': len(conflicts), 'conflicting_samples': len({s for s, p in conflicts}),
        'conflicting_starts': len({p for s, p in conflicts}),
        'target_containing_conflicting_groups_all_traits': sum(any(a['target'] for a in v) for v in conflicts.values()),
        'target_containing_conflicting_groups_nonNC': sum(any(a['target'] and a['trait'] != 'NC' for a in v) for v in conflicts.values()),
        'direct_target_impacts_all_traits': coverage(impacts_all), 'direct_target_impacts_nonNC': coverage(impacts_assay),
        'conflicts_with_single_common_raw_SNV_depth': snv_consistent,
        'conflicts_with_every_MNP_depth_strictly_below_raw_SNV_depth': lower_mnp,
        'conflicts_with_source_lexical_first_allele_raw_SNV': first_snv,
        'normalization_failed_unique_depth_assumption': bool(conflicts),
        'actual_CAO_CAF_WT_q_log_contrasts_or_predictions_decoded': False}
    return summary, table, list(impacts_all.values())


def self_test():
    # Poison both unrequested count/frequency fields. They must not be decoded.
    fields = list(IDENTITY_DP) + ['CAO', 'CAF']
    def stream(dp, ref='A', alt='C', rid='x'):
        vals = dict(record_id=rid, Sample='S', **{'#CHROM': 'toy'}, POS='1', REF=ref, ALT=alt,
                    POS_NT='1', REF_NT='A', ALT_NT='C', DP=str(dp))
        return io.BytesIO(('\t'.join(fields) + '\n').encode() +
            b'\t'.join([vals[k].encode() for k in IDENTITY_DP] + [b'\xffCAO_SEALED', b'\xffCAF_SEALED']) + b'\n')
    contract = {'members': {'S': {'sample': 'S', 'ROI': 1, 'trait': 'R'}},
                'targets': [{'ROI': 1, 'edit_key': '1:A>C', 'hgvs_nt': 'v', 'is_SNV': True}]}
    row1 = list(selected_rows(stream(100), IDENTITY_DP))
    row2 = list(selected_rows(stream(90, 'AA', 'CC', 'y'), IDENTITY_DP))
    status = [{'record_id': 'x', 'VAR_ID': 'toy:1:A:C', 'frequency_pass': 'TRUE', 'globally_rejected': 'FALSE'},
              {'record_id': 'y', 'VAR_ID': 'toy:1:AA:CC', 'frequency_pass': 'TRUE', 'globally_rejected': 'FALSE'}]
    a, _, _ = reconstruct([row1 + row2], [status], contract)
    assert a['conflicting_groups'] == 1 and a['retained_sample_start_groups'] == 1
    assert a['conflicts_with_every_MNP_depth_strictly_below_raw_SNV_depth'] == 1
    assert a['direct_target_impacts_nonNC']['target_IDs'] == 1
    status[1]['globally_rejected'] = 'TRUE'
    b, _, _ = reconstruct([row1 + row2], [status], contract)
    assert b['conflicting_groups'] == 0
    status[1]['globally_rejected'] = 'FALSE'
    row2[0]['DP'] = '100'
    c, _, _ = reconstruct([row1 + row2], [status], contract)
    assert c['conflicting_groups'] == 0
    # Median of an MNP's per-base depths is used before, not after, alleles are combined.
    row2[0]['DP'] = '80'
    extra = {**row2[0], 'record_id': 'z', 'POS_NT': '2', 'DP': '100'}
    status.append({**status[1], 'record_id': 'z'})
    d, _, _ = reconstruct([row1 + row2 + [extra]], [status], contract)
    assert d['counts']['collapsed_alleles_with_per_base_DP_variation'] == 1 and d['conflicting_groups'] == 1
    return {'status': 'PASS_SYNTHETIC_DP_ONLY', 'poisoned_CAO_CAF_unread': True, 'checks': 7}


def main(extract, rf, contract, out):
    out = Path(out)
    if out.exists():
        raise FileExistsError('Fresh failure-only audit directory required')
    extract, rf, contract = Path(extract), Path(rf), Path(contract)
    em, rm = extract / 'sealed_manifest.json', rf / 'sealed_manifest.json'
    if sha(em) != EXTRACT_SEAL or sha(rm) != RF_SEAL or sha(contract) != CONTRACT_SHA:
        raise ValueError('Authorized immutable seal/metadata hash mismatch')
    e, r = [json.loads(x.read_text(encoding='utf8')) for x in (em, rm)]
    if e['phase'] != 'extract' or r['phase'] != 'rf' or r['input_hashes']['extract_manifest'] != EXTRACT_SEAL:
        raise ValueError('Stage identity/lineage mismatch')
    files = {}
    for roi in (1, 2):
        for folder, receipt, name in ((extract, e, 'rows_' + str(roi) + '.tsv'),
                                      (rf, r, 'row_status_' + str(roi) + '.tsv')):
            f = folder / name
            if sha(f) != receipt['files'][name]:
                raise ValueError('Authorized selected input hash mismatch')
            files[str(f)] = sha(f)
    metadata = json.loads(contract.read_text(encoding='utf8'))
    inputs = []
    try:
        for roi in (1, 2):
            inputs.append((extract / ('rows_' + str(roi) + '.tsv')).open('rb'))
        rowiters = [selected_rows(f, IDENTITY_DP) for f in inputs]
        mask_streams = [(rf / ('row_status_' + str(roi) + '.tsv')).open('rb') for roi in (1, 2)]
        inputs.extend(mask_streams)
        maskiters = [selected_rows(f, RF_FIELDS) for f in mask_streams]
        summary, conflicts, impacts = reconstruct(rowiters, maskiters, metadata)
    finally:
        for f in inputs:
            f.close()
    out.mkdir(parents=True)
    for name, table in [('conflicting_groups.csv', conflicts), ('direct_target_impacts.csv', impacts)]:
        with (out / name).open('w', encoding='utf8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(table[0]) if table else ['none'])
            writer.writeheader()
            writer.writerows(table)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf8')
    receipt = {'status': 'COMPLETED_AUTHORIZED_DP_ONLY_FAILURE_AUDIT',
        'code_sha256': sha(__file__), 'selected_fields': list(IDENTITY_DP), 'RF_fields': list(RF_FIELDS),
        'input_hashes': files, 'extract_seal': EXTRACT_SEAL, 'RF_seal': RF_SEAL, 'contract_sha256': CONTRACT_SHA,
        'output_hashes': {p.name: sha(p) for p in out.iterdir()},
        'no_CAO_CAF_normalization_or_downstream_values_decoded': True}
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--extract')
    parser.add_argument('--rf')
    parser.add_argument('--contract')
    parser.add_argument('--out')
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2))
    else:
        if not all((args.extract, args.rf, args.contract, args.out)):
            parser.error('Four explicit selected-input/output paths are required')
        main(args.extract, args.rf, args.contract, args.out)
