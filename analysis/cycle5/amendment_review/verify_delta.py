"""Verify only the four requested pre-outcome clarification paragraphs changed."""
if not __debug__:raise SystemExit('Refuse disabled assertions')
import hashlib
import json
from pathlib import Path

base=Path(__file__).resolve().parent.parent
path=base/'cambray_pilot/amendment_01_draft.md'
new=path.read_bytes();text=new.decode('utf-8');lines=text.splitlines()
old={
31:'3. Omit each of the 16 original related groups in turn, in ascending numeric order. Retain every deletion result, the retained denominators, both-class availability, every baseline comparison and the minimum and maximum across deletions. Recalculate the same metrics without fitting or tuning anything. The full-data criterion remains decisive. These ranges are influence diagnostics, not confidence intervals or performance on new experiments.',
38:'The symmetric additive square with endpoints `(0,0)` and intermediates `(m,-m)` and `(-m,m)` has a joint-only obstruction under every independently bounded per-vertex, per-trait perturbation of magnitude epsilon exactly when `m > 2 epsilon`. The equality case can open a route through weak-threshold ties. More generally, each intermediate-to-endpoint-minimum comparison changes by at most twice the corresponding error bound. This is a deterministic resolution statement for a declared uncertainty set, not a new general theorem or an empirically estimated error bound.',
50:'The primary endpoint is a property of paired measured trait profiles. Independent phenotype shuffling would destroy their covariance. A fixed-endpoint square control that merely exchanges the two intermediate profiles leaves existence of a route unchanged and is not an informative null test. Generic additive antagonism, shared preprocessing, incomplete routes and source uncertainty remain explicit alternatives. No additional correlation model or permuted p-value is introduced without an identifiable design and a separate pre-outcome specification.',
58:'Final outcome categories are: A, a reproducible restriction at the precisely supported scope; B, a descriptive local restriction with stated uncertainty or replication limits; C, a quantitatively supported small-effect exclusion only if calibration and a practical threshold become available; D, an inconclusive test with the demonstrated resolution boundary; or E, a failed source/provenance gate. The current evidence does not authorize C or an independently replicated biological claim under A. Categories D and E are not biological null findings.'}
assert hashlib.sha256(new).hexdigest()=='92c5028dff57df644ad18479f79b4be2e604b7bbe7798e30cdd3dfe5a7fa94d0'
for line,paragraph in old.items():
    assert text.count(lines[line-1])==1
    text=text.replace(lines[line-1],paragraph)
reconstructed=hashlib.sha256(text.encode('utf-8')).hexdigest()
assert reconstructed=='d635d6957550e74bc4551d7f1fc14d7516a67602a58324369984785cd88fbe90',reconstructed
assert hashlib.sha256(path.read_bytes()).hexdigest()==hashlib.sha256(new).hexdigest()
receipt={'status':'passed','original_reviewed_sha256':reconstructed,
    'clarified_amendment_sha256':hashlib.sha256(new).hexdigest(),'changed_paragraph_lines':sorted(old),
    'only_four_requested_paragraphs_changed':True,'biological_sources_opened':0,
    'method':'Restore the four originally reviewed paragraphs in memory and recover the exact original byte hash.',
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
with Path(__file__).with_name('delta_check.json').open('x',encoding='utf-8') as stream:json.dump(receipt,stream,indent=2);stream.write('\n')
print(json.dumps(receipt))
