\name{crit.fun}
\alias{crit.fun}
\title{Find posterior probability threshold to control FDR}
\description{
Find posterior probability threshold to control FDR  
}
\usage{
crit.fun(x, cc) 
}
\arguments{
  \item{x}{ x is one minus the posterior probabilities of being in
    a specific DE pattern. If there is only one DE pattern, then
    x is the posterior probabilities of being EE.}
  \item{cc}{ cc is FDR to be controlled. For example, to control
    FDR at 0.05, set cc=0.05.}
}

\value{
  \code{crit.fun} returns a threshold so that if used in identifying
  genes in a specific DE pattern, FDR can be controlled at cc.
  Those genes with posterior probability of being in that
  specific DE pattern greater than this threshold are claimed to be
  in that specific DE pattern.
}
\author{Ming Yuan, Ping Wang, Deepayan sarkar, Michael Newton, and Christina Kendziorski}
\references{
  Newton, M.A., Noueiry, A., Sarkar, D., and Ahlquist, P. (2004).
  Detecting differential gene expression with a semiparametric
  hierarchical mixture model. \emph{Biostatistics} \bold{5}, 155-176.
}
\examples{
data(gould)
pattern <- ebPatterns(c("1,1,1,0,0,0,0,0,0,0",
                        "1,2,2,0,0,0,0,0,0,0"))
gg.em.out <- emfit(gould, family = "GG", hypotheses = pattern, num.iter = 10)
gg.post.out <- postprob(gg.em.out, gould)$pattern 
gg.crit <- crit.fun(gg.post.out[,1],0.05)
# number of DE genes 
sum(gg.post.out[,2] > gg.crit)

pattern4 <- ebPatterns(c("1, 1, 1, 1, 1, 1, 1, 1, 1, 1", 
                         "1, 2, 2, 2, 2, 2, 2, 2, 2, 2", 
                         "1,2,2,1,1,1,1,1,2,2", 
                         "1,1,1,1,1,1,1,1,2,2"))
gg4.em.out <- emfit(gould, family = "GG", pattern4, num.iter = 10)
gg4.post.out <- postprob(gg4.em.out, gould)$pattern
gg4.crit <- crit.fun(1-gg4.post.out[,2], 0.05)
# number of genes in pattern 2, a DE pattern
sum(gg4.post.out[,2] > gg4.crit)
}
\keyword{}
