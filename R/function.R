# x: posterior probability of EE; cc: FDR cutoff
# genes with posterior probability of DE greater than
# the returned value is called DE at FDR level cc. 
crit.fun <- function(x,cc){
  yy <- cumsum(sort(x))/(1:length(x))
  mm <- yy<cc
  index <- sum(mm)
  if(index>0){
    out <- 1-sort(x)[index]}
  if(index==0){
    out <- 1
  }
  names(out) <- NULL
  return(out)
}

