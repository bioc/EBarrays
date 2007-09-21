\name{ebarraysFamily-class}
\alias{ebarraysFamily-class}
\alias{eb.createFamilyGG}
\alias{eb.createFamilyLNN}
\alias{eb.createFamilyLNNMV}
\alias{coerce,character,ebarraysFamily-method}
\alias{show,ebarraysFamily-method}
\title{Class of Families to be used in the EBarrays package}
\description{
  Objects used as family in the \code{emfit} function. 

  The package contains three functions that create such objects for the
  three most commonly used families, Gamma-Gamma, Lognormal-Normal and
  Lognormal-Normal with modified variances. Users may create their own 
  families as well.
}
\usage{
eb.createFamilyGG()
eb.createFamilyLNN()
eb.createFamilyLNNMV()
}
\details{
  The \code{\link{emfit}} function can potentially fit models
  corresponding to several different Bayesian conjugate families. This
  is specified as the \code{family} argument, which ultimately has to be
  an object of formal class ``ebarraysFamily'' with some specific slots
  that determine the behavior of the `family'.

  For users who are content to use the predefined GG, LNN and LNNMV models, no
  further details than that given in the documentation for
  \code{\link{emfit}} are necessary. If you wish to create your own
  families, read on.
}
\value{
  Objects of class ``ebarraysFamily'' for the three predefined families
  Gamma-Gamma , Lognormal-Normal and Lognormal-Normal with modified 
  variances. 
}
\section{Objects from the Class}{
  Objects of class ``ebarraysFamily'' can be created by calls of the
  form \code{new("ebarraysFamily", ...)}. Predefined objects
  corresponding to the GG, LNN and LNNMV models can be created by
  \code{eb.createFamilyGG()} , \code{eb.createFamilyLNN()} and
  \code{eb.createFamilyLNNMV()}. The same
  effect is achieved by coercing from the strings \code{"GG"}, \code{"LNN"}
  and \code{"LNNMV"} by \code{as("GG", "ebarraysFamily")}, \code{as("LNN",
  "ebarraysFamily")} and \code{as("LNNMV", "ebarraysFamily")}.
}
\section{Slots}{
  An object of class ``ebarraysFamily'' extends the class
  \code{"character"} (representing a short hand name for the class) and
  should have the following slots (for more details see the source
  code):

  \describe{
    \item{\code{description}:}{
      A not too long character string describing the family
    }
    \item{\code{link}:}{
      function that maps user-visible parameters to the parametrization that
      would be used in the optimization step (e.g. \code{log(sigma^2)}
      for LNN). This allows the user to think in terms of familiar
      parametrization that may not necessarily be the best when
      optimizing w.r.t. those parameters.
    }
    \item{\code{invlink}:}{
      inverse of the link function
    }
    \item{\code{thetaInit}:}{
      function of a single argument \code{data} (matrix containing raw
      expression values), that calculates and returns as a numeric
      vector initial estimates of the parameters (in the parametrization
      used for optimization)
    }
    \item{\code{f0}:}{
      function taking arguments \code{theta} and a list called
      \code{args}. \code{f0} calculates the negative log likelihood at
      the given parameter value \code{theta} (again, in the
      parametrization used for optimization). This is called from
      \code{emfit}. When called, only genes with positive intensities
      across all samples are used.
    }
    \item{\code{f0.pp}:}{
      \code{f0.pp} is essentially the same as \code{f0} except the terms
      common to the numerator and denominator when calculating posterior
      odds may be removed. It is called from \code{postprob}.
    }
    \item{\code{f0.arglist}:}{
      function that takes arguments \code{data}, \code{patterns} (of
      class ``ebarraysPatterns'') and \code{groupid} (for LNNMV family
      only) and returns a list with two components, \code{common.args} and
      \code{pattern.args}. \code{common.args} is a list of arguments to
      \code{f0} that don't change from one pattern to another, whereas
      \code{pattern.args[[i]][[j]]} is a similar list of arguments, but
      specific to the columns in \code{pattern[[i]][[j]]}. Eventually,
      the two components will be combined for each pattern and used as
      the \code{args} argument to \code{f0}.
    }
    \item{\code{logDensity}:}{
      function of two arguments \code{x} (data vector, containing log
      expressions) and \code{theta} (parameters in user-visible
      parametrization). Returns log marginal density of the natural log
      of intensity for the corresponding theoretical model. Used in
      \code{plotMarginal}
    }
    \item{\code{lower.bound}:}{
      vector of lower bounds for the argument \code{theta} of
      \code{f0}. Used in \code{optim}
    }
    \item{\code{upper.bound}:}{
      vector of upper bounds for the argument \code{theta} of
      \code{f0}. 
    }
  }
}
\author{Ming Yuan, Ping Wang, Deepayan Sarkar, Michael Newton, and Christina Kendziorski}
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
  \code{\link{emfit}}, \code{\link{optim}}, \code{\link{plotMarginal}}
}
\examples{
show(eb.createFamilyGG())
show(eb.createFamilyLNN())
show(eb.createFamilyLNNMV())
}
\keyword{models}



