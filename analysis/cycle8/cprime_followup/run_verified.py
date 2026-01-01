"""Safe entry point retaining the exact reviewed analysis source."""
if not __debug__:
    raise RuntimeError('Optimized Python is prohibited before any data access')
from pathlib import Path
import hashlib
import argparse
import importlib.util
import json

p=Path(__file__).with_name('analysis.py')
if hashlib.sha256(p.read_bytes()).hexdigest()!='45f75f86ef97325ae17bcabd306504dc57f5d885e81d3244a71eb43a5870619f':
    raise RuntimeError('The reviewed analysis source changed')
parser=argparse.ArgumentParser(description='Exact reviewed calculation in a fresh directory; outputs are under OUT/results_v1')
parser.add_argument('--input',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
if args.out.exists():raise FileExistsError('Refusing existing output')
if hashlib.sha256(args.input.read_bytes()).hexdigest()!='eeeb40c167e689ae4266bbe49fdef7c8711163ca426b217adaee14f0b4a6afdf':
    raise ValueError('Changed B1 source')
protocol=p.with_name('protocol.md').read_bytes()
if hashlib.sha256(protocol).hexdigest()!='a3acdaad3beedba8793352fe402a587850c66a92bfb6c33b7ee42ca539c1decd':
    raise ValueError('Changed exploratory protocol')
spec=importlib.util.spec_from_file_location('reviewed_cprime_followup',p)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
# Only historical filesystem constants change. All scientific functions and
# the actual source/protocol bytes remain pinned to the reviewed run.
module.SOURCE=args.input.resolve()
module.HERE=args.out.resolve()
args.out.mkdir(parents=True)
(args.out/'protocol.md').write_bytes(protocol)
module.main()
(args.out/'wrapper_receipt.json').write_text(json.dumps({'wrapper_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_sha256':module.SOURCE_SHA,'filesystem_overrides_only':['SOURCE','HERE']},indent=2)+'\n',encoding='utf-8')
