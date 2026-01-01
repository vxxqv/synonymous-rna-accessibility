"""Declare every within-site pair from the already authenticated metadata only."""
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
CONTRACT = ROOT / 'work/research_cycle8/comt_pilot/implementation_v3/selective_contract.json'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    d = json.loads(CONTRACT.read_text(encoding='utf-8'))
    rows = sorted(d['targets'], key=lambda x: x['hgvs_nt'].encode('utf-8'))
    assert len(rows) == len({r['hgvs_nt'] for r in rows}) == 125
    grouped = defaultdict(list)
    for r in rows:
        grouped[int(r['native_codon_site'])].append(r)
    pairs = []
    for site, rr in sorted(grouped.items()):
        for a, b in itertools.combinations(rr, 2):
            ea = {int(x.split(':')[0]): x.split('>')[1] for x in a['edit_key'].split(';')}
            eb = {int(x.split(':')[0]): x.split('>')[1] for x in b['edit_key'].split(';')}
            reference = {}
            for r in (a, b):
                for x in r['edit_key'].split(';'):
                    p, change = x.split(':'); base = change.split('>')[0]
                    assert int(p) not in reference or reference[int(p)] == base
                    reference[int(p)] = base
            hamming = sum(ea.get(p, reference[p]) != eb.get(p, reference[p]) for p in reference)
            pairs.append({'a': a['hgvs_nt'], 'b': b['hgvs_nt'], 'site': site, 'ROI': a['ROI'],
                          'both_SNV': a['is_SNV'] and b['is_SNV'],
                          'observed_edit_string_distance': hamming,
                          'same_edit_start': min(ea) == min(eb),
                          'same_start_is_not_authenticated_shared_normalization': True})
    summary = {}
    for label, vv, pp in [
        ('all', rows, pairs),
        ('SNV', [r for r in rows if r['is_SNV']], [p for p in pairs if p['both_SNV']])
    ]:
        per = Counter(int(r['native_codon_site']) for r in vv)
        summary[label] = {'variants': len(vv), 'sites': len(per), 'pairs': len(pp),
                          'sites_with_pairs': sum(n >= 2 for n in per.values()),
                          'singleton_sites': sum(n == 1 for n in per.values()),
                          'pair_hamming_histogram': dict(sorted(Counter(p['observed_edit_string_distance'] for p in pp).items())),
                          'same_edit_start_pairs': sum(p['same_edit_start'] for p in pp),
                          'site_variant_count_histogram': dict(sorted(Counter(per.values()).items()))}
    output = {'status': 'METADATA_ONLY_CANDIDATE_UNIVERSE', 'input_sha256': sha(CONTRACT),
              'source_sha256': sha(Path(__file__)), 'summaries': summary,
              'sites': [{'site': s, 'ROI': rr[0]['ROI'], 'ids': [r['hgvs_nt'] for r in rr]} for s, rr in sorted(grouped.items())],
              'pairs': pairs, 'outcome_values_read': False,
              'interpretation': 'These are complete pair comparisons within the selected allele set, not complete sequence space or evolutionary routes.'}
    target = OUT / 'metadata_inventory.json'
    if target.exists(): raise FileExistsError(target)
    target.write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'summary': summary, 'receipt_sha256': sha(target)}))

if __name__ == '__main__':
    if not __debug__: raise SystemExit('Optimized Python unsupported')
    main()
