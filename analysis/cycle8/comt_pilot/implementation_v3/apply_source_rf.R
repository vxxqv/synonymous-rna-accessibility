# Executable only through the guarded Python driver. No network or package install.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop('Expected pilot_root prepared_dir output_dir')
pilot <- normalizePath(args[1],winslash='/',mustWork=TRUE)
prepared <- normalizePath(args[2],winslash='/',mustWork=TRUE)
out <- args[3]
if(dir.exists(out)) stop('Fresh RF output required')
source(file.path(pilot,'source_helpers.R'))
h <- load_comt_helpers(pilot);e<-h$env
dir.create(out,recursive=TRUE)
features<-c('VAR_ID','UP_REF_NT','DOWN_REF_NT','TYPE','log10_CAF',
            'R2_MINUS_MED_RP','R1_PLUS_MED_BQ','R2_MINUS_MED_BQ','R2_MINUS_MED_NM')
all_rejects<-data.table::data.table(ROI=integer(),VAR_ID=character(),false_gdna_samples=integer())
metadata_names<-c('VAR_ID','TYPE','Sample')
for(roi in 1:2){
  sim<-data.table::fread(file.path(prepared,paste0('simulation_',roi,'.tsv')),sep='\t',check.names=FALSE)
  truth<-data.table::fread(file.path(prepared,paste0('truth_',roi,'.tsv')),sep='\t',check.names=FALSE)
  sim[,TYPE:=e$get_var_type(REF,ALT),by=seq_len(nrow(sim))]
  if(roi==1L) sim<-sim[TYPE=='SNP' | log10(CAF)>-4]
  train<-e$prep_train_dt(sim,features,truth_dt=truth)
  expected<-sim[LOCATION=='CDS' & CAF<.3]
  if(nrow(train)!=nrow(expected))stop('Simulation matrix row loss')
  if(nrow(train)<10L || length(unique(train$Truth))!=2L)stop('Source ten-fold RF design unavailable')
  if(any(!is.finite(as.matrix(train[,-'Truth']))))stop('Nonfinite RF training matrix')
  folds<-e$create_kfolds(train,kfolds=10)
  tuning<-train[!folds[[1]]]
  cv<-e$run_nested_cv(tuning,model='rf',hyperparam_prop=.3,kfolds=5)
  mtry<-round(median(cv[[2]]$mtry))
  if(!is.finite(mtry)||!mtry%in%3:10)stop('Unexpected source mtry')
  set.seed(9)
  fitted<-randomForest::randomForest(Truth~.,train,mtry=mtry)
  rows<-data.table::fread(file.path(prepared,paste0('rows_',roi,'.tsv')),sep='\t',check.names=FALSE)
  prediction<-e$prep_train_dt(rows,features)
  eligible<-rows[LOCATION=='CDS' & CAF<.3]
  eligible_ids<-paste(eligible[['#CHROM']],eligible$POS,eligible$REF,eligible$ALT,sep=':')
  if(nrow(prediction)!=nrow(eligible) || !identical(as.character(prediction$VAR_ID),eligible_ids))stop('Prediction identity or row loss')
  wanted<-setdiff(names(train),'Truth')
  actual<-setdiff(names(prediction),metadata_names)
  if(!identical(actual,wanted))stop('Training/prediction feature mismatch')
  if(any(!is.finite(as.matrix(prediction[,..actual]))))stop('Nonfinite prediction matrix')
  probabilities<-predict(fitted,prediction,type='prob')
  if(!identical(colnames(probabilities),c('FALSE','TRUE')))stop('Unexpected class labels')
  prediction[,Prediction:=factor(probabilities[,2]>=.5,levels=c('FALSE','TRUE'))]
  prediction[,probability_TRUE:=probabilities[,2]]
  prediction[,record_id:=eligible$record_id]
  gdna<-paste('COMT',if(roi==1L)c(37:39,112:115) else 116:119,'R',sep='_')
  fp<-prediction[Sample%in%gdna & Prediction=='FALSE',.N,by=.(Sample,VAR_ID,TYPE)]
  counts<-fp[TYPE=='SNP' & N==1,.N,by=VAR_ID]
  rejected<-counts[N>=if(roi==1L)4 else 2,VAR_ID]
  selected<-counts[VAR_ID%in%rejected,.(ROI=roi,VAR_ID,false_gdna_samples=N)]
  all_rejects<-rbind(all_rejects,selected)
  rows[,VAR_ID:=paste(`#CHROM`,POS,REF,ALT,sep=':')]
  rows[,frequency_pass:=log10(CAF)>-5.8]
  rows[,globally_rejected:=VAR_ID%in%rejected]
  data.table::fwrite(rows[,.(record_id,VAR_ID,frequency_pass,globally_rejected)],file.path(out,paste0('row_status_',roi,'.tsv')),sep='\t')
  data.table::fwrite(prediction[,.(record_id,Sample,VAR_ID,TYPE,Prediction,probability_TRUE)],file.path(out,paste0('predictions_',roi,'.tsv')),sep='\t')
  saveRDS(list(fitted_model=fitted,tuning=cv,initial_folds=folds,mtry=mtry,
      source_seed=9,feature_names=wanted,source_training_rows=nrow(train),
      interpretation='New fixed-runtime preprocessing; not historical state'),
      file.path(out,paste0('model_',roi,'.rds')),version=3)
}
data.table::fwrite(all_rejects,file.path(out,'rejections.tsv'),sep='\t')
writeLines('PASS_NEW_FIXED_RUNTIME_RF_APPLICATION',file.path(out,'status.txt'))
