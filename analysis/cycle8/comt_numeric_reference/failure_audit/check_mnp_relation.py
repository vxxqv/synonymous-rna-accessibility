"""Read only the DP-only conflict CSV produced by the frozen probe."""
if not __debug__:
    raise RuntimeError('Optimized Python prohibited')
import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def main(source, out):
    source, out = Path(source), Path(out)
    if out.exists():
        raise FileExistsError(out)
    counts = {'groups': 0, 'all_MNP_at_or_below_common_SNV': 0,
              'at_least_one_MNP_strictly_below_common_SNV': 0,
              'some_MNP_equal_common_SNV': 0, 'any_MNP_above_common_SNV': 0}
    with source.open(newline='', encoding='utf8') as f:
        for row in csv.DictReader(f):
            counts['groups'] += 1
            alleles = json.loads(row['all_depth_and_identity_records'])
            s = {F(a['DP']) for a in alleles if a['raw_SNV']}
            m = {F(a['DP']) for a in alleles if not a['raw_SNV']}
            if len(s) != 1 or not m:
                raise AssertionError('Common-SNV comparison not defined')
            d = next(iter(s))
            counts['all_MNP_at_or_below_common_SNV'] += all(v <= d for v in m)
            counts['at_least_one_MNP_strictly_below_common_SNV'] += any(v < d for v in m)
            counts['some_MNP_equal_common_SNV'] += d in m
            counts['any_MNP_above_common_SNV'] += any(v > d for v in m)
    result = {'counts': counts, 'scope': 'Independent DP-only CSV relation check; no additional source fields read',
              'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    main(a.input, a.out)
