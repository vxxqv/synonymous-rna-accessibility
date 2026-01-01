"""Format-only adapter for pinned deposited files; frozen scientific core unchanged."""
from pathlib import Path
import hashlib
import json
import sys
import count_analysis as core

ORIGINAL=core.parse_counts
HEADER=b'Genotype\tOccurrences\tMutationCount\n'

def parse_headerless(raw,expected_git_blob=None):
    if expected_git_blob:
        got=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        if got!=expected_git_blob:raise ValueError('Original source blob mismatch')
    first=raw.decode('utf-8-sig').splitlines()[0].split()
    if len(first)!=3 or first[0]!=core.WT:
        raise ValueError('Amendment applies only to verified headerless files starting with native WT')
    return ORIGINAL(HEADER+raw,None)

def main():
    here=Path(__file__).resolve().parent
    amendment=json.loads((here/'format_freeze.json').read_text())
    for name,sha in amendment['files'].items():
        if hashlib.sha256((here/name).read_bytes()).hexdigest()!=sha:raise RuntimeError('Format amendment changed: '+name)
    core.parse_counts=parse_headerless
    core.main()

if __name__=='__main__':main()
