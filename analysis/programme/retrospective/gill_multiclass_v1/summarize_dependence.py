"""Summarize already frozen normalization comparisons without fitting or new outcomes."""
from pathlib import Path
import argparse
import csv
import hashlib
import json

SOURCE_SHA256 = 'cdbdd83d6861c54d1cd5094990cbeb306803349c63bdd0519caf391676c819f6'

def summarize(repo, out):
    source = repo / 'results/programme/gill_multiclass_v1/validate/validation.json'
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA256:
        raise RuntimeError('Frozen validation source differs')
    if out.exists():
        raise RuntimeError('Refusing existing output')
    data = json.loads(raw)
    order = data['scenario_order']
    pairs = data['decomposition_pair_order']
    if len(order) != 9 or len(pairs) != 45:
        raise RuntimeError('Unexpected frozen comparison inventory')
    rows = []
    results = []
    for pseudocount, partitions in data['all_summaries'].items():
        for partition, summary in partitions.items():
            buckets = {'different_RNA_shared_DNA': [], 'different_RNA_different_DNA': [],
                       'same_RNA_different_DNA': [], 'identical_scenario': []}
            values = summary['scenario_residual_pooled_cosines']
            if len(values) != 45:
                raise RuntimeError('Incomplete comparison list')
            for (a, b), value in zip(pairs, values, strict=True):
                ra, da = order[a]
                rb, db = order[b]
                if a == b:
                    category = 'identical_scenario'
                elif ra == rb:
                    category = 'same_RNA_different_DNA'
                elif da == db:
                    category = 'different_RNA_shared_DNA'
                else:
                    category = 'different_RNA_different_DNA'
                buckets[category].append(value)
                rows.append({'partition': partition, 'pseudocount': pseudocount,
                             'RNA_a': ra + 1, 'DNA_a': da + 8,
                             'RNA_b': rb + 1, 'DNA_b': db + 8,
                             'category': category, 'pooled_residual_cosine': value})
            if {k: len(v) for k, v in buckets.items()} != {
                'different_RNA_shared_DNA': 9, 'different_RNA_different_DNA': 18,
                'same_RNA_different_DNA': 9, 'identical_scenario': 9}:
                raise RuntimeError('Comparison partition is incomplete')
            results.append({'partition': partition, 'pseudocount': pseudocount,
                            'groups': summary['class_count'],
                            'categories': {k: {'n': len(v), 'minimum': min(v), 'maximum': max(v)}
                                           for k, v in buckets.items()},
                            'raw_rank_agreement': summary['rank_agreement_summary']})
    out.mkdir(parents=True, exist_ok=False)
    with (out / 'all_normalization_comparisons.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {'status': 'descriptive_summary_of_frozen_comparisons',
              'source_sha256': SOURCE_SHA256, 'comparisons': len(rows), 'settings': results,
              'interpretation': 'Comparisons share libraries and genotypes. They are dependent diagnostics, not biological replicates or calibrated uncertainty.',
              'selection': 'Every original partition, pseudocount and pair retained. Category ranges summarized after validation; not a new prospective endpoint.',
              'anchor_invariance': 'A scalar reference subtraction within each class is annihilated by Q because Q1=0. It also preserves within-class ranks and endpoint-relative barriers. No shared reference anchor was used by the frozen analysis.'}
    (out / 'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({'status': report['status'], 'comparisons': len(rows),
                      'primary_validation': next(x for x in results if x['partition'] == 'validation' and x['pseudocount'] == '0.5')}, indent=2))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    summarize(args.repo_root.resolve(), args.out.resolve())
