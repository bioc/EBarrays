\name{ebplots}
\alias{ebplots}
\alias{checkCCV}
\alias{checkModel}
\alias{checkVarsQQ}
\alias{checkVarsMar}
\alias{plotMarginal}
\alias{plotCluster}
\alias{plot.ebarraysEMfit}
\title{Various plotting routines in the EBarrays package}
\description{
  Various plotting routines, used for diagnostic purposes
}
\usage{
checkCCV(data, useRank = FALSE, f = 1/2)
checkModel(data, fit, model = c("gamma", "lognormal", "lnnmv"),
           number = 9, nb = 10, cluster = 1, groupid = NULL)
checkVarsQQ(data, groupid, \dots)
checkVarsMar(data, groupid, xlab, ylab, \dots)
plotMarginal(fit, data, kernel = "rect", n = 100,
             bw = "nrd0", adjust = 1, xlab, ylab,\dots)
plotCluster(fit, data, cond = NULL, ncolors = 123, sep=TRUE,
            transform=NULL)

\S3method{plot}{ebarraysEMfit}(x, data, plottype="cluster", \dots)
}
\arguments{
  \item{data}{ data, as a ``matrix'' or ``ExpressionSet''}
  \item{useRank}{ logical. If \code{TRUE}, ranks of means and c.v.-s are
    used in the scatterplot
  }
  \item{f}{ passed on to \code{\link{lowess}}}
  \item{fit, x}{ object of class ``ebarraysEMfit'', typically produced by a
    call to \code{\link{emfit}}
  }
  \item{model}{ which theoretical model use for Q-Q plot. Partial string
    matching is allowed
  }
  \item{number}{ number of bins for checking model assumption.}
  \item{nb}{ number of data rows included in each bin for checking model
    assumption}
  \item{cluster}{ check model assumption for data in that cluster }
  \item{groupid}{ an integer vector indicating which group each sample
    belongs to. groupid for samples not included in the analysis should
    be 0. }
  \item{kernel, n, bw, adjust}{ passed on to \code{\link{density}}}
  \item{cond}{ a vector specifying the condition for each replicate}
  \item{ncolors}{ different number of colors in the plot}
  \item{xlab, ylab}{ labels for x-axis and y-axis}
  \item{sep}{ whether or not to draw horizontal lines between clusters}
  \item{transform}{ a function to transform the original data in
    plotting}
  \item{plottype}{ a character string specifying the type of the plot.
    Available options are "cluster" and "marginal". The default plottype
    "cluster" employs function 'plotCluster' whereas the "marginal"
    plottype uses function 'plotMarginal'.}
  \item{\dots}{ extra arguments are passed to the \code{qqmath},
    \code{histogram} and \code{xyplot} call used to produce the final result
  }
}
\details{
  \code{checkCCV} checks the constant coefficient of variation assumption
  made in the GG and LNN models.
  \code{checkModel} generates QQ plots for subsets of (log) intensities
  in a small window. They are used to check the Log-Normal assumption on
  observation component of the LNN and LNNMV models and the Gamma
  assumption on observation component of the GG model.
  \code{checkVarsQQ} generates QQ plot for gene specific sample
  variances. It is used to check the assumption of a scaled inverse
  chi-square prior on gene specific variances, made in the LNNMV model.
  \code{checkVarsMar} is another diagnostic tool to check this
  assumption. The density histogram of gene specific sample variances
  and the density of the scaled inverse chi-square distribution with
  parameters estimated from data will be plotted.
  \code{checkMarginal} generates predictive marginal distribution from
  fitted model and compares with estimated marginal (kernel) density of
  data. Available for the GG and LNN models only. 
  \code{plotCluster} generate heatmap for gene expression data with clusters
}
\value{
  \code{checkModel}, \code{checkVarsQQ} and \code{checkVarsMar}
  return an object of class
  ``trellis'', using function in the Lattice package. Note that in
  certain situations, these may need to be explicitly `print'-ed to have
  any effect.
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
\seealso{\code{\link{emfit}}, \code{\link{lowess}}}
\examples{
}
\keyword{models}
