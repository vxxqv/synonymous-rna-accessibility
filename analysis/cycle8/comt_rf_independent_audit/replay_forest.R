# Independent feature reconstruction and saved-state prediction. No training.
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4L)stop('Expected pilot_root prepared_dir rf_dir NEW_output_dir')
pilot<-normalizePath(a[1],winslash='/',mustWork=TRUE)
prepared<-normalizePath(a[2],winslash='/',mustWork=TRUE)
rf<-normalizePath(a[3],winslash='/',mustWork=TRUE)
out<-a[4]
if(dir.exists(out))stop('New output required')
.libPaths(c(file.path(pilot,'runtime','library'),.Library))
suppressPackageStartupMessages(library(randomForest))
if(as.character(getRversion())!='4.6.1' || as.character(packageVersion('randomForest'))!='4.7.1.2')stop('Runtime mismatch')
if(Sys.setlocale('LC_COLLATE','C')!='C')stop('C collation unavailable')
options(contrasts=c('contr.treatment','contr.poly'),digits=17)
dir.create(out,recursive=TRUE)
readtab<-function(path)read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,colClasses='character',na.strings=character(),quote='',comment.char='')
type<-function(ref,alt){
  if(nchar(ref)!=nchar(alt))stop('Only substitution source records expected')
  if(nchar(ref)==1L)return('SNP')
  if(nchar(ref)==2L)return('di_nt_MNP')
  if(nchar(ref)==3L){
    hd<-sum(strsplit(ref,'',fixed=TRUE)[[1]]!=strsplit(alt,'',fixed=TRUE)[[1]])
    if(hd==2L && !grepl(',',ref,fixed=TRUE))return('di_nt_MNP')
    if(hd==2L && grepl(',',ref,fixed=TRUE))return('2_nt_HAPLO')
    return('tri_nt_MNP')
  }
  # Source fixed TYPE factor does not admit longer haplotypes. Never recode.
  return('outside_source_TYPE_levels')
}
reports<-list()
for(roi in 1:2){
  rows<-readtab(file.path(prepared,paste0('rows_',roi,'.tsv')))
  if(anyDuplicated(rows$record_id))stop('Duplicate extracted identity')
  freq<-as.numeric(rows$CAF)
  if(any(!is.finite(freq)) || any(freq<0|freq>1))stop('Invalid CAF')
  use<-rows$LOCATION=='CDS' & freq<.3
  variant<-paste(rows[['#CHROM']],rows$POS,rows$REF,rows$ALT,sep=':')
  types<-mapply(type,rows$REF,rows$ALT,USE.NAMES=FALSE)
  if(any(!types[use]%in%c('SNP','di_nt_MNP','tri_nt_MNP')))stop('Unsupported source TYPE')
  # Source factor levels are discovered on all rows before CAF restriction.
  # Encode explicitly instead of importing prep_train_dt/model.matrix.
  mat<-matrix(numeric(),nrow=sum(use),ncol=0)
  for(name in c('UP_REF_NT','DOWN_REF_NT')){
    levels<-sort(unique(rows[[name]]),method='radix')
    if(length(levels)<2L)stop('Source factor contrasts need two levels')
    for(level in levels[-1L]){
      mat<-cbind(mat,as.numeric(rows[[name]][use]==level))
      colnames(mat)[ncol(mat)]<-paste0(name,level)
    }
  }
  for(level in c('di_nt_MNP','tri_nt_MNP')){
    mat<-cbind(mat,as.numeric(types[use]==level));colnames(mat)[ncol(mat)]<-paste0('TYPE',level)
  }
  numeric_names<-c('log10_CAF','R2_MINUS_MED_RP','R1_PLUS_MED_BQ','R2_MINUS_MED_BQ','R2_MINUS_MED_NM')
  for(name in numeric_names){
    values<-if(name=='log10_CAF')log10(freq[use]) else as.numeric(rows[[name]][use])
    if(any(!is.finite(values)))stop('Nonfinite model feature')
    mat<-cbind(mat,values);colnames(mat)[ncol(mat)]<-name
  }
  saved<-readRDS(file.path(rf,paste0('model_',roi,'.rds')))
  if(!identical(colnames(mat),saved$feature_names))stop('Independent feature columns disagree')
  if(saved$source_seed!=9 || saved$fitted_model$ntree!=500L)stop('Saved seed/tree metadata mismatch')
  if(saved$mtry!=round(median(saved$tuning[[2]]$mtry)))stop('Saved mtry aggregation mismatch')
  probs<-predict(saved$fitted_model,newdata=as.data.frame(mat,check.names=FALSE),type='prob')
  if(!identical(colnames(probs),c('FALSE','TRUE')))stop('Unexpected forest classes')
  actual<-readtab(file.path(rf,paste0('predictions_',roi,'.tsv')))
  ids<-rows$record_id[use]
  if(!identical(ids,actual$record_id) || !identical(variant[use],actual$VAR_ID) || !identical(rows$Sample[use],actual$Sample) || !identical(types[use],actual$TYPE))stop('Full prediction identity/order mismatch')
  observed<-as.numeric(actual$probability_TRUE)
  difference<-probs[,2]-observed
  label<-ifelse(probs[,2]>=.5,'TRUE','FALSE')
  trace<-data.frame(record_id=ids,VAR_ID=variant[use],Sample=rows$Sample[use],TYPE=types[use],replayed_probability=probs[,2],saved_probability=observed,difference=difference,replayed_class=label,saved_class=actual$Prediction,check.names=FALSE)
  write.table(trace,file.path(out,paste0('trace_',roi,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  if(any(!is.finite(observed)) || any(abs(difference)>1e-14) || any(label!=actual$Prediction))stop('Probability or class mismatch; trace retained')
  gdna<-paste0('COMT_',if(roi==1L)c(37:39,112:115) else 116:119,'_R')
  failed<-data.frame(Sample=rows$Sample[use],VAR_ID=variant[use],TYPE=types[use],Prediction=label,stringsAsFactors=FALSE)
  failed<-failed[failed$Sample%in%gdna & failed$Prediction=='FALSE',]
  counts<-if(nrow(failed))aggregate(rep(1L,nrow(failed)),by=list(Sample=failed$Sample,VAR_ID=failed$VAR_ID,TYPE=failed$TYPE),FUN=sum) else data.frame(Sample=character(),VAR_ID=character(),TYPE=character(),x=integer())
  names(counts)[4]<-'n'
  counts<-counts[counts$TYPE=='SNP' & counts$n==1L,]
  totals<-table(counts$VAR_ID)
  rejected<-names(totals)[totals>=if(roi==1L)4L else 2L]
  mask<-readtab(file.path(rf,paste0('row_status_',roi,'.tsv')))
  if(!identical(mask$record_id,rows$record_id) || !identical(mask$VAR_ID,variant))stop('Full status identity/order mismatch')
  fpass<-log10(freq)>-5.8
  gpass<-variant%in%rejected
  if(!identical(mask$frequency_pass,ifelse(fpass,'TRUE','FALSE')) || !identical(mask$globally_rejected,ifelse(gpass,'TRUE','FALSE')))stop('Full row status mismatch')
  rej<-readtab(file.path(rf,'rejections.tsv'));rej<-rej[as.integer(rej$ROI)==roi,]
  if(anyDuplicated(rej$VAR_ID) || !setequal(rej$VAR_ID,rejected))stop('Rejection identity mismatch')
  if(length(rejected) && any(as.integer(rej$false_gdna_samples)!=as.integer(totals[rej$VAR_ID])))stop('Rejection count mismatch')
  reports[[roi]]<-data.frame(ROI=roi,extracted_rows=nrow(rows),predicted_rows=nrow(trace),rejected_SNPs=length(rejected),frequency_pass=sum(fpass),max_probability_difference=max(abs(difference)),class_mismatches=sum(label!=actual$Prediction),saved_training_rows=saved$source_training_rows,source_seed_metadata=saved$source_seed,ntree=saved$fitted_model$ntree,mtry=saved$mtry)
}
write.table(do.call(rbind,reports),file.path(out,'summary.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines('PASS_SAVED_FOREST_TRACE_NOT_TRAINING_VALIDATION',file.path(out,'status.txt'))
