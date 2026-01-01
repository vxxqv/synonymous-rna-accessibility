"""Create a path/decompression-only adapter; preserve the original verifier."""
from pathlib import Path
import ast
import difflib
import hashlib
import json

BASE=Path(__file__).resolve().parent
def sha(raw):return hashlib.sha256(raw).hexdigest()

def main():
    source=BASE/'verify_result.py';raw=source.read_bytes()
    if sha(raw)!='e971da5fab874f0843b58d2e85dc9fc6b28939254eb4ee3c19d73466f46bc46f':raise RuntimeError('Original verifier drift')
    text=raw.decode('utf-8');adapted=text
    changes=[
      ("parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path,required=True)",
       "parser=argparse.ArgumentParser();parser.add_argument('--repo','--root',dest='repo',type=Path,required=True)"),
      ("parser.add_argument('--result',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()",
       "parser.add_argument('--result',type=Path);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()\n    result_path=args.result if args.result is not None else args.repo/'results/cycle7/robustness.json.gz'"),
      ("input_hashes={str(frozen_path):FREEZE_SHA,str(args.result):RESULT_SHA}",
       "input_hashes={str(frozen_path):FREEZE_SHA}"),
      ("result_raw=read_pinned(args.result,RESULT_SHA)",
       "saved_result=result_path.read_bytes()\n    result_raw=gzip.decompress(saved_result) if result_path.suffix=='.gz' else saved_result\n    need(sha(result_raw)==RESULT_SHA,'Decoded successor result hash mismatch')\n    input_hashes[str(result_path)]=sha(saved_result)"),
    ]
    for old,new in changes:
        if adapted.count(old)!=1:raise RuntimeError('Expected adapter target not unique')
        adapted=adapted.replace(old,new)
    old_functions={n.name:ast.dump(n,include_attributes=False) for n in ast.parse(text).body if isinstance(n,ast.FunctionDef)}
    new_functions={n.name:ast.dump(n,include_attributes=False) for n in ast.parse(adapted).body if isinstance(n,ast.FunctionDef)}
    changed=[name for name in old_functions if old_functions[name]!=new_functions[name]]
    if changed!=['main']:raise RuntimeError('Unexpected algorithm modification')
    path=BASE/'verify_portable.py'
    if path.exists():raise RuntimeError('Refusing existing portable verifier')
    path.write_bytes(adapted.encode('utf-8'))
    (BASE/'portable.diff').write_text(''.join(difflib.unified_diff(text.splitlines(True),adapted.splitlines(True),fromfile='verify_result.py',tofile='verify_portable.py')),encoding='utf-8')
    receipt={'status':'PORTABLE_ADAPTER_READY','original_sha256':sha(raw),'portable_sha256':sha(path.read_bytes()),
             'unchanged_function_asts':[name for name in old_functions if name not in changed],
             'changed_functions':changed,
             'scope':'Only CLI root alias, optional default successor path, gzip decoding and file-level hash bookkeeping changed.',
             'decoded_successor_sha256':'861b7061db34532698b074c22409aceb16925a56348b5dcdbf313bcb2dbd40aa'}
    (BASE/'portable_preparation.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))

if __name__=='__main__':main()
