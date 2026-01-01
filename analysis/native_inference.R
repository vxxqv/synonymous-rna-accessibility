args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "."
d <- read.csv(file.path(root,"results","native_paired.csv"),na.strings="")
deduplicated <- length(args)>1 && args[2]=="deduplicate"
if(deduplicated) d <- d[!duplicated(d$locus_tag),,drop=FALSE]
set.seed(1729)
groups <- split(seq_len(nrow(d)),d$transcript_id)
cluster_boot <- function(response, adjust=FALSE) {
  use <- !is.na(d[[response]])
  dd <- d[use,,drop=FALSE]
  if(nrow(dd)==0) return(c(estimate=NA,lower=NA,upper=NA,n=0,clusters=0))
  gg <- split(seq_len(nrow(dd)),dd$transcript_id)
  x <- if(adjust) model.matrix(~d_pairs+d_span+d_internal,dd) else matrix(1,nrow(dd),1)
  y <- dd[[response]]
  cross <- lapply(gg,function(i) list(xx=crossprod(x[i,,drop=FALSE]),xy=crossprod(x[i,,drop=FALSE],y[i])))
  fit <- function(weights) {
    xx <- matrix(0,ncol(x),ncol(x)); xy <- matrix(0,ncol(x),1)
    for(j in which(weights>0)) {xx <- xx+weights[j]*cross[[j]]$xx;xy <- xy+weights[j]*cross[[j]]$xy}
    tryCatch(as.numeric(solve(xx,xy)[1]),error=function(e) NA_real_)
  }
  estimate <- fit(rep(1,length(gg)))
  b <- replicate(10000,fit(tabulate(sample.int(length(gg),length(gg),replace=TRUE),nbins=length(gg))))
  ci <- quantile(b,c(.025,.975),na.rm=TRUE,names=FALSE)
  c(estimate=estimate,lower=ci[1],upper=ci[2],n=nrow(dd),clusters=length(gg))
}
results <- rbind(log_effective_raw=cluster_boot("d_log_effective"),
 log_effective_adjusted=cluster_boot("d_log_effective",TRUE),
 log_abundance_raw=cluster_boot("d_log_abundance"),
 excess_paired_nonempty=cluster_boot("d_excess"),
 native_pair_retention=cluster_boot("d_retention"))
write.csv(results,file.path(root,"results",if(deduplicated) "native_inference_unique_loci.csv" else "native_inference.csv"),row.names=TRUE)
capture.output(sessionInfo(),file=file.path(root,"results","R_session.txt"))
print(results)
