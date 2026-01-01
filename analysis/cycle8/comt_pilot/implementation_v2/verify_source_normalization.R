# Entirely invented fixture, using the unchanged author collapse helper.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Expected pilot_root fresh_output')
if(file.exists(args[2]))stop('Fresh output required')
source(file.path(args[1],'source_helpers.R'))
e<-load_comt_helpers(args[1])$env
d<-data.table::data.table(Sample=factor(rep('S',3)),POS=c(1L,1L,1L),
  REF=c('AA','AA','A'),ALT=c('CC','CC','G'),POS_NT=c(1L,2L,1L),
  REF_NT=c('A','A','A'),ALT_NT=c('C','C','G'),FILTER='PASS',
  UP_REF_NT='A',DOWN_REF_NT='C',AA_POS='1',REF_AA='Q',ALT_AA='P',
  LOCATION='CDS',DP=c(100L,102L,101L),CAO=c(10L,12L,30L),CAF=c(.1,.12,.3))
d[,`#CHROM`:='toy']
c<-e$postprocess_vcf_summary_dt(d)
stopifnot(nrow(c)==2L,length(unique(c$DP))==1L)
c[,WT:=unique(DP)-sum(CAO)]
c[,q:=(CAO+.5)/(WT+.5)]
data.table::fwrite(c[,.(VAR_ID,DP,CAO,WT,q)],args[2],sep='\t')
cat('PASS_INVENTED_SOURCE_COLLAPSE\n')
