\name{utilities}
\alias{utilities}
\alias{ebPatterns}
\alias{show,ebarraysPatterns-method}
\alias{ebarraysPatterns-class}
\title{Utility functions for the EBarrays package}
\description{
  Utilitiy functions for the EBarrays package
}
\usage{
ebPatterns(x, ordered=FALSE)
}
\arguments{
  \item{x}{
    x can be a character vector (of length > 2) (see example), or an
    arbitrary connection which should provide patterns, one line for
    each pattern. If \code{x} is a character vector of length 1, it is
    assumed to be the name of a file (since there's no point in a
    patterns object with only one pattern) which is then opened and
    treated as a connection.
  }
  \item{ordered}{ logical variable specifying whether the pattern is
    ordered or not}
}
\details{
  \code{ebPatterns} creates objects that represent a collection of
  hypotheses to be used by \code{emfit}.
}
\value{
  \code{ebPatterns} creates an Object of class ``ebarraysPatterns'', to
  be used in other functions such as \code{\link{emfit}}. This is
  nothing more than a list (and can be treated as such as far as
  indexing goes) and is used only for method dispatch.
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
\seealso{\code{\link{emfit}}}
\examples{
patterns <- ebPatterns(c("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1",
                         "1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2
                          2 2 2"), TRUE)
show(patterns)
}
\keyword{models}
