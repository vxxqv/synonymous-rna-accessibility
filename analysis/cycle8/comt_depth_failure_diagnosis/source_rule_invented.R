# Invented values only. Unchanged collapse helper plus literal notebook WT code.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Expected pilot root and fresh output directory')
if(dir.exists(args[2]))stop('Fresh output required')
source(file.path(args[1],'source_helpers.R'))
e<-load_comt_helpers(args[1])$env
dir.create(args[2],recursive=TRUE)
d<-data.table::data.table(Sample=factor(rep('S',4)),ROI=factor('ROI_1'),
 Input=factor('RNA'),Fraction=factor('total'),Quadrant=factor('none'),POS=1L,
 REF=c('A','A','AA','AA'),ALT=c('C','G','CC','CC'),POS_NT=c(1L,1L,1L,2L),
 REF_NT='A',ALT_NT=c('C','G','C','C'),FILTER='PASS',UP_REF_NT='A',DOWN_REF_NT='C',
 AA_POS='1',REF_AA='Q',ALT_AA='P',LOCATION='CDS',
 DP=c(100L,100L,90L,90L),CAO=c(10L,20L,5L,5L),CAF=c(.1,.2,5/90,5/90))
d[,`#CHROM`:='toy']
keys<-c('Sample','ROI','Input','Fraction','Quadrant','VAR_ID')
c<-e$postprocess_vcf_summary_dt(d,key_features=keys)
c2<-e$postprocess_vcf_summary_dt(d[4:1],key_features=keys)
stopifnot(isTRUE(all.equal(c[,.(VAR_ID,DP,CAO)],c2[,.(VAR_ID,DP,CAO)])),nrow(c)==3L)
assign_source<-function(x){
 w<-x[,.(CAO_sum=sum(CAO),DP=unique(DP)),by=.(ROI,Input,Fraction,Quadrant,Sample,POS)]
 w[,WT_CAO:=DP-CAO_sum]
 x[,WT_key:=paste(Sample,POS,sep=':')];w[,WT_key:=paste(Sample,POS,sep=':')]
 x[,WT_CAO:=w[match(x[,WT_key],WT_key),WT_CAO]]
 list(assigned=x,denominators=w)
}
r<-assign_source(data.table::copy(c))
stopifnot(nrow(r$denominators)==2L,identical(r$denominators$DP,c(100,90)),
 all(r$assigned$WT_CAO==65),identical(r$denominators$WT_CAO,c(65,55)))
reversed<-assign_source(data.table::copy(c)[3:1])
stopifnot(all(reversed$assigned$WT_CAO==55))
data.table::fwrite(r$assigned[,.(VAR_ID,DP,CAO,WT_CAO)],file.path(args[2],'assigned.tsv'),sep='\t')
data.table::fwrite(r$denominators,file.path(args[2],'multiple_denominators.tsv'),sep='\t')
writeLines(c('PASS_INVENTED_SOURCE_RULE','No unique post-collapse depth requirement exists in the executed source expression.',
 'Two distinct DP values emit two WT rows; match chooses the first.',
 'Source helper key sorting makes the supplied raw-row reversal immaterial in this fixture.',
 'Manually reversing collapsed rows changes first-match WT from65to55.'),file.path(args[2],'receipt.txt'))
cat('PASS_INVENTED_SOURCE_RULE\n')
