from pathlib import Path
import hashlib,json,csv,os
P=Path(__file__).resolve().parent
H=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
files={}
root=P/'runtime/library'; ext=Path('\\\\?\\'+str(root.resolve()))
for base,dirs,names in os.walk(ext):
 for name in names:
  p=Path(base)/name;files[str(p.relative_to(ext)).replace('\\','/')]=H(p)
versions=list(csv.DictReader((P/'runtime/installed_versions.tsv').open(),delimiter='\t'))
manifest={'status':'INSTALLED_NEW_RUNTIME_SYNTHETIC_SMOKE_PASSED_NOT_HISTORICAL','R_version':'4.6.1','Rscript':'C:/Program Files/R/R-4.6.1/bin/Rscript.exe','Rscript_sha256':H(Path('C:/Program Files/R/R-4.6.1/bin/Rscript.exe')),'RNGkind':['Mersenne-Twister','Inversion','Rejection'],'source_set_seed':9,'LC_COLLATE':'C','data_table_threads':1,'foreach_backend':'sequential','caret_default_resampling':{'method':'boot','number':25,'selectionFunction':'best','allowParallel_flag':True,'actual_workers':1},'RF_ntree':500,'versions':versions,'package_files':files,'download_manifest_sha256':H(P/'runtime/download_manifest.json'),'source_helpers_sha256':H(P/'source_helpers.R'),'helper_source_sha256':{n:H(P/n) for n in ['modeling_utils.r','summarization_utils.r']},'helper_AST_receipt_sha256':H(P/'runtime/helper_ast_hashes.tsv'),'synthetic_RF_receipt_sha256':H(P/'runtime/synthetic_rf_receipt.rds'),'source_count_rows_read':0,'binary_authentication_limit':'Downloaded over official CRAN HTTPS; version/identity and internal MD5 files verified; archive SHA256 recorded. No separately published cryptographic binary digest was available in the binary package index.','install_notes':['Initial R installer failed on a long Windows path; its local partial directory and log are preserved.','Same downloaded binaries extracted with extended paths into task-local runtime/library; packages loaded and synthetic RF repeated exactly.','No global or system packages modified.']}
(P/'runtime/runtime_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Pinned',len(files),'package files; runtime SHA',H(P/'runtime/runtime_manifest.json'))
