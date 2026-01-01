"""Fetch only the known public calibration source, mask every numeric token.

Never execute author code or request counts, PSI tables or outcome missingness.
The unmasked byte response is retained as sealed source evidence and not printed.
"""
from pathlib import Path
import hashlib
import json
import re
import requests

HERE = Path(__file__).resolve().parent
URL = 'https://raw.githubusercontent.com/lehner-lab/Scaling_Law/b87064c39a553b0759d86c2b428813c58dc17ae2/003_HEK293_Phylogeny_Library_Analysis/002_psi_scores.R'
EXPECTED_GIT_BLOB = 'c3ee712c600e926a4a653410f860ac49c10b08bd'
response = requests.get(URL, timeout=40)
receipt = {'url': URL, 'status': response.status_code,
           'sha256': hashlib.sha256(response.content).hexdigest(),
           'bytes': len(response.content),
           'scope': 'author code only; numerical literals masked before display'}
if response.status_code != 200:
    (HERE/'retrieval_failure.json').write_text(json.dumps(receipt,indent=2),encoding='utf8')
    raise RuntimeError('Source access failed; no retry')
blob = hashlib.sha1(b'blob '+str(len(response.content)).encode()+b'\0'+response.content).hexdigest()
receipt.update(git_blob=blob,expected_git_blob=EXPECTED_GIT_BLOB)
if blob != EXPECTED_GIT_BLOB:
    raise RuntimeError('Source differs from known pinned author blob')
(HERE/'calibration_source_sealed.R').write_bytes(response.content)
text = response.content.decode('utf8')
# Conservative: also mask digits embedded in quoted strings, comments and names.
# This obscures some indices but prevents numeric outcome literals reaching display.
masked, n = re.subn(r'\d+(?:\.\d*)?(?:[eE][+-]?\d+)?', '<NUM>', text)
(HERE/'calibration_source_masked.txt').write_text(masked,encoding='utf8')
receipt['masked_numeric_tokens'] = n
receipt['masked_sha256'] = hashlib.sha256(masked.encode()).hexdigest()
(HERE/'calibration_source_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
print(json.dumps(receipt,indent=2))
