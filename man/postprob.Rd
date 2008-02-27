\name{postprob}
\alias{postprob}
\alias{postprob,ebarraysEMfit,ExpressionSet-method}
\alias{postprob,ebarraysEMfit,matrix-method}
\alias{show,ebarraysPostProb-method}
\alias{ebarraysPostProb-class}
\title{Calculates posterior probabilities for expression patterns}
\description{
  Takes the output from emfit and calculates the posterior probability
  of each of the hypotheses, for each gene.
}
\usage{
postprob(fit, data, ...)
}
\arguments{
  \item{fit}{ output from \code{\link{emfit}}}
  \item{data}{ a numeric matrix or an object of class ``ExpressionSet''
    containing the data, typically the same one used in the \code{emfit}
    fit supplied below.
  }
  \item{\dots}{ other arguments, ignored}
}
\value{
  An object of class ``ebarraysPostProb''. Slot \code{joint} is an three
  dimensional array of probabilities. Each element gives the posterior
  probability that a gene belongs to certain cluster and have certain
  pattern. \code{cluster} is a matrix of probabilities with  number of
  rows given by the number of genes in \code{data} and as many
  columns as the number of clusters for the fit. \code{pattern} is a
  matrix of probabilities with  number of rows given by the number of
  genes in \code{data} and as many columns as the number of patterns for
  the fit. It additionally contains a slot `hypotheses' containing
  these hypotheses.
}
\author{Ming Yuan, Ping Wang, Deepayan Sarkar, Michael Newton, and
  Christina Kendziorski}
\references{
  
  Newton, M.A., Kendziorski, C.M., Richmond, C.S., Blattner, F.R. (2001).
  On differential variability of expression ratios: Improving
  statistical inference about gene expression changes from microarray data.
  Journal of Computational Biology 8:37-52.

  Kendziorski, C.M., Newton, M.A., Lan, H., Gould, M.N. (2003).
  On parametric empirical Bayes methods for comparing multiple groups
  using replicated gene expression profiles.
  Statistics in Medicine 22:3899-3914.

  Newton, M.A. and Kendziorski, C.M.
  Parametric Empirical Bayes Methods for Microarrays in
  The analysis of gene expression data: methods and software. Eds.
  G. Parmigiani, E.S. Garrett, R. Irizarry and S.L. Zeger,
  New York: Springer Verlag, 2003.

  Newton, M.A., Noueiry, A., Sarkar, D., and Ahlquist, P. (2004).
  Detecting differential gene expression with a semiparametric
  hierarchical mixture model. Biostatistics 5: 155-176.

  Yuan, M. and Kendziorski, C. (2006). A unified approach for simultaneous
  gene clustering and differential expression identification.
  Biometrics 62(4): 1089-1098.

}
\seealso{
  \code{\link{emfit}}
}
\examples{
data(sample.ExpressionSet) ## from Biobase
eset <- exprs(sample.ExpressionSet)
patterns <- ebPatterns(c("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1",
                         "1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2"))
gg.fit <- emfit(data = eset, family = "GG", hypotheses = patterns, verbose = TRUE)
prob <- postprob(gg.fit,eset)
}
\keyword{models}
