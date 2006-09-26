\name{ebplots}
\alias{ebplots}
\alias{checkCCV}
\alias{checkConvergence}
\alias{plotMarginal}
\alias{checkModel}
\title{Various plotting routines in the EBarrays package}
\description{
  Various plotting routines, used for diagnostic purposes
}
\synopsis{
plotMarginal(fit, data, kernel = "rect", n = 100, bw = "nrd0", adjust = 1,
             xlab, ylab, \dots)
}
\usage{
checkCCV(data, useRank = FALSE, f = 1/2)
plotMarginal(fit, data, kernel = "rect", n = 100,
             bw = "nrd0", adjust = 1, \dots)
checkModel(data, model = c("gamma", "lognormal"),
           number = 9, nb = 10)
checkConvergence(\dots, dropfirst = FALSE)
}
\arguments{
  \item{data}{ data, as a \code{"matrix"}, \code{"exprSet"} or \code{"ExpressionSet"}}
  \item{useRank}{ logical. If \code{TRUE}, ranks of means and c.v.-s are
    used in the scatterplot
  }
  \item{f} {passed to \code{\link{lowess}}}
  \item{fit}{ object of class ``ebarraysEmFit'', typically produced by a
    call to \code{\link{emfit}}
  }
  \item{kernel, n, bw, adjust}{ passed on to \code{\link{density}}}
  \item{model}{ which theoretical model use for Q-Q plot. Partial string
    matching is allowed
  }
  \item{number}{ number of subsets or bins (panels in the Trellis
    display) for which Q-Q plot is drawn
  }
  \item{nb}{ number of observations per bin (panel) }
  \item{\dots}{
    In the case of \code{checkConvergence}, matrices or vectors (treated
    as column vectors) with same number of rows. These are combined into
    a single matrix, then each column is scaled between 0 and 1 and
    plotted simultaneously.

    For \code{plotMarginal}, extra arguments are passed to the
    \code{xyplot} call used to produce the final result.
  }
  \item{dropfirst}{ logical. If \code{TRUE}, first row will be dropped }
}
\details{

  Needs to be written
}
\value{
  \code{plotMarginal} and \code{checkModel} return an object of class
  ``trellis'', using function in the Lattice package. Note that in
  certain situations, these may need to be explicitly `print'-ed to have
  any effect.
}
\seealso{\code{\link{emfit}}, \code{\link{lowess}}}
\author{Christina Kendziorski, Michael Newton and Deepayan Sarkar}
\examples{
data(gould)
checkCCV(gould)

}
\keyword{models}



