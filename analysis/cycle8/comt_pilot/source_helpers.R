# New fixed-runtime adapter. This file never opens assay or simulation data.
load_comt_helpers <- function(base_dir) {
  base_dir <- normalizePath(base_dir, winslash='/', mustWork=TRUE)
  .libPaths(c(file.path(base_dir,'runtime','library'), .Library))
  for(p in c('data.table','caret','randomForest','e1071'))
    library(p,character.only=TRUE)
  if(as.character(getRversion())!='4.6.1') stop('R version mismatch')
  versions <- c(data.table='1.18.6.1',caret='7.0.1',randomForest='4.7.1.2',e1071='1.7.17')
  for(p in names(versions))
    if(as.character(packageVersion(p))!=versions[[p]]) stop('Package version mismatch')
  source_hashes <- c(
    'modeling_utils.r'='219e33c36a0f03dca82b6d504819091f8e36c015489e60544cb901b9bd741a1d',
    'summarization_utils.r'='f0a741b087e9e06e6d32cd947f2faf2dcedd59272d3e0a87921de118bad740bd')
  for(f in names(source_hashes))
    if(digest::digest(file=file.path(base_dir,f),algo='sha256')!=source_hashes[[f]])
      stop('Helper source hash mismatch')
  foreach::registerDoSEQ()
  data.table::setDTthreads(1L)
  if(foreach::getDoParWorkers()!=1L) stop('Sequential runtime required')
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  if(Sys.setlocale('LC_COLLATE','C')!='C') stop('Cannot set fixed collation')
  options(contrasts=c('contr.treatment','contr.poly'))
  wanted <- list(
    'summarization_utils.r'=c('get_var_type','annotate_calls','numeric_median',
      'collapse_mnps','summarize_str_list','postprocess_vcf_summary_dt'),
    'modeling_utils.r'=c('cont_table_factor_levels','get_accuracy',
      'get_precision','get_sensitivity','get_recall','create_kfolds',
      'prep_train_dt','hyperparam_tune','run_nested_cv'))
  env <- new.env(parent=globalenv())
  evidence <- list()
  for(f in names(wanted)) {
    expr <- parse(file.path(base_dir,f), keep.source=FALSE)
    names_found <- character()
    for(e in expr) {
      if(!is.call(e) || !as.character(e[[1]]) %in% c('<-','=')) next
      name <- as.character(e[[2]])
      if(length(name)!=1L || !name %in% wanted[[f]]) next
      if(name %in% names_found) stop('Duplicate requested helper')
      eval(e,envir=env)
      if(is.call(e[[3]]) && identical(e[[3]][[1]],as.name('function'))) {
        if(!identical(formals(env[[name]]),e[[3]][[2]]) ||
           !identical(body(env[[name]]),e[[3]][[3]])) stop('Changed helper AST')
      }
      names_found <- c(names_found,name)
      evidence[[paste(f,name,sep=':')]] <- digest::digest(e,algo='sha256')
    }
    if(!setequal(names_found,wanted[[f]])) stop('Missing requested helper')
  }
  if(!identical(env$cont_table_factor_levels,c('FALSE','TRUE'))) stop('Class order')
  # Non-RF model branches remain in unchanged function bodies but are forbidden
  # by the caller. No unrelated package library() expression is evaluated.
  list(env=env, ast_hashes=evidence)
}
