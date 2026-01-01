"""Synthetic format-adapter audit. No deposited count file is opened."""
from pathlib import Path
import hashlib
import importlib.util
import json
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[2]
REPO=BASE/'outputs/synonymous-rna-accessibility'
SOURCE=REPO/'analysis/cycle6/fas_pilot'
EXPECTED_ADAPTER='11addfe6ee9c02f7fc17b9b6a414167cc0d98baccbab71bf7fa899fd7b713bb6'
FROZEN_COMMIT='7d8faf59fddde21ef4f1177f13ab5d2144d192e6'


def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
def require(x,msg):
 if not x:raise RuntimeError(msg)


def main():
 output=HERE/'receipt.json';require(not output.exists(),'No overwrite')
 adapter_raw=(SOURCE/'run_headerless.py').read_bytes()
 require(sha(adapter_raw)==EXPECTED_ADAPTER,'Adapter drift')
 freeze=json.loads((SOURCE/'freeze.json').read_text())
 frozen_checks=[]
 for name,expected in freeze['frozen_files'].items():
  path=SOURCE/name;raw=path.read_bytes()
  old=subprocess.check_output(['git','-C',str(REPO),'show',f'{FROZEN_COMMIT}:analysis/cycle6/fas_pilot/{name}'])
  require(raw==old and sha(raw)==expected,'Frozen core changed: '+name)
  frozen_checks.append({'name':name,'sha256':expected,'bytes_identical_to_original_commit':True})
 original_freeze=subprocess.check_output(['git','-C',str(REPO),'show',f'{FROZEN_COMMIT}:analysis/cycle6/fas_pilot/freeze.json'])
 require(original_freeze==(SOURCE/'freeze.json').read_bytes(),'Original freeze mutated')
 cfg=json.loads((SOURCE/'source_contract.json').read_text())
 require(len(cfg['files'])==12 and all(re.fullmatch('[0-9a-f]{40}',d['git_blob']) for d in cfg['files']),
         'CLI source contract does not pin twelve source blobs')
 sys.path.insert(0,str(SOURCE))
 spec=importlib.util.spec_from_file_location('adapter_under_review',SOURCE/'run_headerless.py')
 adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
 core=adapter.core;tests=[]
 wt=core.WT;others=[s for s in cfg['sequences'] if s!=wt][:2]
 def row(s,c,m=None):
  if m is None:m=sum(a!=b for a,b in zip(s,core.ANCESTOR))
  return f'{s}\t{c}\t{m}'
 lines=[row(wt,17),row(others[0],0),row(others[1],10**45)]
 cases={'LF':('\n'.join(lines)+'\n').encode(),
        'CRLF':('\r\n'.join(lines)+'\r\n').encode(),
        'no_trailing_newline':'\n'.join(lines).encode(),
        'blank_inner_row':(lines[0]+'\n\n'+lines[1]+'\n'+lines[2]+'\n').encode(),
        'spaces':('\n'.join(s.replace('\t','   ') for s in lines)+'\n').encode()}
 expected={wt:17,others[0]:0,others[1]:10**45}
 for name,raw in cases.items():
  require(adapter.parse_headerless(raw,blob(raw))==expected,'First/all rows not preserved: '+name)
  require(adapter.parse_headerless(raw,blob(raw))==core.parse_counts(adapter.HEADER+raw,None),'Core parser equality failed')
  tests.append({'name':name,'status':'same_three_records_including_first_WT_and_zero'})
 # Instrument only the original parser reference, without running analysis.
 calls=[];saved=adapter.ORIGINAL
 def capture(raw,expected_git_blob):calls.append((raw,expected_git_blob));return expected
 adapter.ORIGINAL=capture
 raw=cases['CRLF'];require(adapter.parse_headerless(raw,blob(raw))==expected,'Spy result')
 require(calls==[(adapter.HEADER+raw,None)],'Adapter changed original row bytes or parser interface')
 adapter.ORIGINAL=saved
 tests.append({'name':'exact_byte_forwarding','status':'Only fixed header prepended; every original byte retained'})
 def reject(name,raw,expected_hash=None,substring=None):
  try:adapter.parse_headerless(raw,expected_hash)
  except (ValueError,IndexError,UnicodeError) as error:
   if substring:require(substring in str(error),'Unexpected refusal reason: '+name)
   tests.append({'name':name,'status':'refused','exception':type(error).__name__,'message':str(error)})
  else:raise RuntimeError('Bad source accepted: '+name)
 reject('bad_git_blob_before_utf8_decode',b'\xff','0'*40,'Original source blob mismatch')
 reject('modified_source_with_preedit_hash',cases['LF'].replace(b'\t17\t',b'\t18\t'),blob(cases['LF']),'Original source blob mismatch')
 reject('duplicate_native_row',(lines[0]+'\n'+lines[0]+'\n').encode())
 reject('duplicate_later_variant',(lines[0]+'\n'+lines[1]+'\n'+lines[1]+'\n').encode())
 reject('negative_count',(row(wt,-1)+'\n').encode())
 reject('fractional_count',(row(wt,'1.5')+'\n').encode())
 reject('unicode_count',(row(wt,'\u0661')+'\n').encode())
 reject('wrong_ancestor_annotation',(row(wt,1,999)+'\n').encode())
 reject('wrong_first_identity',(row(others[0],1)+'\n'+lines[0]+'\n').encode())
 reject('malformed_later_row',(lines[0]+'\n'+others[0]+'\t2\n').encode())
 reject('outside_design_later_row',(lines[0]+'\n'+row('N'+wt[1:],2)+'\n').encode())
 reject('all_zero_library',(row(wt,0)+'\n'+row(others[0],0)+'\n').encode())
 reject('already_headered',adapter.HEADER+cases['LF'])
 reject('empty_source',b'')
 # BOM is not stripped from the original row after prefixing; refusal is fail-closed.
 reject('unexpected_BOM',b'\xef\xbb\xbf'+cases['LF'])
 # Calling with a required source hash must never bypass the frozen parser.
 require(adapter.ORIGINAL is core.parse_counts,'Adapter mutated core at import time')
 require(sha((SOURCE/'run_headerless.py').read_bytes())==EXPECTED_ADAPTER,'Source changed during audit')
 result={'status':'PASS_SYNTHETIC_FORMAT_ADAPTER','adapter_sha256':EXPECTED_ADAPTER,
         'amendment_sha256':sha((SOURCE/'format_amendment_v1.md').read_bytes()),
         'original_frozen_commit':FROZEN_COMMIT,'frozen_file_checks':frozen_checks,
         'original_freeze_byte_identical':True,'tests':tests,'test_count':len(tests),
         'real_count_files_opened':False,'biological_outcomes_evaluated':False,
         'CLI_execution_pending':'format_freeze.json not yet present at inspection; root must freeze adapter/amendment and invoke unchanged core guards',
         'qualifications':['Headerless contract requires first native WT record; empty/BOM/wrong-first files refuse.',
                           'Direct parser expected hash is optional for synthetic use; frozen CLI supplies all twelve pinned hashes.',
                           'Secondary source formula is corrected-headerless reconstruction, not authenticated historical result.'],
         'review_script_sha256':sha(Path(__file__).read_bytes())}
 output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'status':result['status'],'tests':len(tests),'frozen_files':len(frozen_checks),'receipt_sha256':sha(output.read_bytes())}))


if __name__=='__main__':main()
