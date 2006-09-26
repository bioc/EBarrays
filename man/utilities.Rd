\name{utilities}
\alias{utilities}
\alias{createExprSet}
\alias{ebPatterns}
\alias{show,ebarraysPatterns-method}
\alias{ebarraysPatterns-class}
\title{Utility functions for the EBarrays package}
\description{
  Utilitiy functions for the EBarrays package
}
\usage{
createExprSet(datafile, repfile)
ebPatterns(x)
}
\arguments{
  \item{datafile}{ file with data}
  \item{repfile}{ file with replication information}
  \item{x}{
    x can be a character vector (of length > 2) (see example), or an
    arbitrary connection which should provide patterns, one line for
    each pattern. If \code{x} is a character vector of length 1, it is
    assumed to be the name of a file (since there's no point in a
    patterns object with only one pattern) which is then opened and
    treated as a connection.
  }
}
\details{
  \code{ebPatterns} creates objects that represent a collection of
  hypotheses to be used by \code{emfit}.

  It, along with \code{createExprSet}, also provides back compatibility
  with an earlier version of \code{EBarrays}, which required input data
  to be in a very specific format. For example, the (earlier) usage
  \code{cleanup(NULL, "data.txt", "repfile.txt", "patternfile.txt")}
  can be replaced by \code{createExprSet("data.txt", "repfile.txt")} and
  \code{ebPatterns("patternfile.txt")}.

  Unless one has data in this particular format, the
  \code{createExprSet} function can be safely ignored.
}

\value{
  \code{createExprSet} creates an object of class \code{"exprSet"}.

  \code{ebPatterns} creates an Object of class ``ebarraysPatterns'', to
  be used in other functions such as \code{\link{emfit}}. This is
  nothing more than a list (and can be treated as such as far as
  indexing goes) and is used only for method dispatch.

}
\seealso{\code{\link{emfit}}}
\author{Christina Kendziorski, Michael Newton and Deepayan Sarkar}
\examples{
patterns <- ebPatterns(c("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1",
                         "1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2"))
show(patterns)
}
\keyword{models}



