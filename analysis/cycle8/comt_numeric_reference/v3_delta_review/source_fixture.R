# Invented rows only. Original helper bodies and literal original WT/match algebra.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L || file.exists(args[3])) stop('pilot_root invented_csv fresh_output required')
Sys.setlocale('LC_COLLATE','C')
source(file.path(args[1],'source_helpers.R'))
e <- load_comt_helpers(args[1])$env
input <- data.table::fread(args[2])
all <- list()
for (fixture in unique(input$fixture)) {
  index <- which(input[['fixture']]==fixture)
  d <- data.table::copy(input[index,])
  d[,fixture:=NULL]
  for(n in c('Sample','ROI','Input','Fraction','Quadrant')) d[,(n):=as.factor(get(n))]
  collapsed <- e$postprocess_vcf_summary_dt(d,
    key_features=c('Sample','ROI','Input','Fraction','Quadrant','VAR_ID'))
  denom <- collapsed[,.(CAO_sum=sum(CAO),DP=unique(DP)),
    by=.(ROI,Input,Fraction,Quadrant,Sample,POS)]
  denom[,WT_CAO:=DP-CAO_sum]
  collapsed[,WT_key:=paste(Sample,POS,sep=':')]
  denom[,WT_key:=paste(Sample,POS,sep=':')]
  collapsed[,WT_CAO:=denom[match(collapsed[,WT_key],WT_key),WT_CAO]]
  collapsed[,q:=(CAO+0.5)/(WT_CAO+0.5)]
  out <- collapsed[,.(fixture=fixture,Sample,POS,VAR_ID,DP,CAO,WT_CAO,q)]
  all[[fixture]] <- out
}
data.table::fwrite(data.table::rbindlist(all),args[3],sep='\t')
cat('PASS_ORIGINAL_HELPER_AND_LITERAL_SOURCE_RULE_INVENTED_ONLY\n')
