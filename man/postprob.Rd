\name{postprob}
\alias{postprob}
\alias{postprob,ebarraysEmFit,exprSet-method}
\alias{postprob,ebarraysEmFit,matrix-method}
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
  \item{data}{ a nunmeric matrix or an object of class ``exprSet''
    containing the data, typically the same one used in the \code{emfit}
    fit supplied below.
  }
  \item{\dots}{ other arguments, ignored}
}
\value{
  An object of class ``ebarraysPostProb'', which is essentially a matrix
  of probabilities with number of rows given by the number of genes in
  \code{data} and as many columns as the number of hypotheses for the
  fit. It additionally contains a slot `hypotheses' containing
  these hypotheses.
}
\seealso{
  \code{\link{emfit}}}
}
\author{ Christina Kendziorski, Michael Newton and Deepayan Sarkar}
\keyword{models}
