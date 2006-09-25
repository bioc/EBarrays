
### Copyright 2003-2004, Christina Kendziorski
### <kendzior@biostat.wisc.edu>, Michael Newton
### <newton@biostat.wisc.edu> and Deepayan Sarkar
### <deepayan@stat.wisc.edu>
###
### This file is part of the EBarrays library for R.  It is made
### available under the terms of the GNU General Public License,
### version 2, or at your option, any later version, incorporated
### herein by reference.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details. 
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
### 02111-1307, USA


checkCCV <- function(data, useRank = FALSE, f = 1/2)
{
    cv <- function(x) {
        return(sqrt(var(x))/mean(x))
    }
    if (is(data, "exprSet") || is(data, "ExpressionSet"))
        data <- exprs(data)
    means <- apply(data, 1, mean)
    cvs <- apply(data, 1, cv)

    if (useRank)
    {
        plot(rank(means), rank(cvs), pch=".",
             xlab=c("Rank of Mean"),
             ylab=c("Rank of CV") )
        lines(lowess(rank(means), rank(cvs), f = f), col = "red")
    }
    else
    {
        plot(means, cvs, pch = ".",
             xlab = "Mean",
             ylab = "CV",
             log="xy")
        lines(lowess(means, cvs, f = f), col="red")
    }
    invisible()
}






checkConvergence <- function(..., dropfirst = FALSE)
    ## x is a matrix. Each column will be scaled between 0 and 1 and
    ## plotted simultaneously
{
    x <- as.matrix(cbind(...))
    if (ncol(x) < 1) return()

    for (i in 1:ncol(x))
    {
        rng <- if (dropfirst) range(x[-1, i])
        x[,i] <- (x[,i] - rng[1]) / diff(rng)
    }
    matplot(t(x))
}









plotMarginal <-
function (fit, data, kernel = "rect", n = 100, bw = "nrd0", adjust = 1,
          xlab, ylab,
          ...)
{

    ## Input: data is a matrix or exprSet/ExpressionSet. fit is the
    ## output from emfit

    ## Output: Plot comparing the marginal distribution under the
    ## theoreetical model to the empirical distribution, on the log
    ## scale.

    ## Note: x_ij values have the same marginal, but x_i1, ... x_im
    ## are not necessarily independent. Still OK to estimate density
    ## from all data, since the result can be thought of as the
    ## average of several density estimates, one for each column.

    if (!require(lattice, quietly = TRUE)) stop("The lattice package could not be loaded")

    if (is(data, "exprSet") || is(data, "ExpressionSet")) data <- exprs(data)
    data <- log(data[data > 0])

    lims <- range(data)

    ## By default, plot average shifted histogram of log expressions
    ## instead of plain histogram. Other kernels also possible

    empmarg <- density(data, bw = "nrd0", adjust = 1,
                       kernel = kernel, n = n,
                       from = lims[1], to = lims[2])

    supp <- empmarg$x
    theta <- fit@family@invlink(fit@thetaEst)
    logmarg <- fit@family@logDensity(supp, theta)

    if (missing(xlab)) xlab <- "log(expressions)"
    if (missing(ylab)) ylab <- "Density"

    xyplot(exp(logmarg) + empmarg$y ~ empmarg$x, type = 'l',
           allow.mult = TRUE,
           key = simpleKey(text = c("Marginal Density of log Expressions from Fitted Model",
                           "Empirical Kernel Density of log Expressions"),
           lines = TRUE, points = FALSE),
           xlab = xlab,
           ylab = ylab, ...)
}











checkModel <-
    function(data, model = c("gamma", "lognormal"),
             number = 9, nb = 10)
{
    if (!require(lattice, quietly = TRUE)) stop("The lattice package could not be loaded")
    model <- match.arg(model)
    if (is(data, "exprSet") || is(data, "ExpressionSet")) data <- exprs(data)
    
    means <- apply(data, 1, mean)
    means <- means[means > 0]
    mean.ranks <- equal.count(rank(means), number = number,
                              overlap = (number - length(means) / nb) / (number - 1))

    if (model == "lognormal")
    {
        ans <- 
            qqmath(~ log(means) | mean.ranks, distribution = qnorm, 
                   scales = list(relation = "free", draw = FALSE),
                   panel = function(x, distribution, ...) {
                       panel.qqmathline(x, distribution = distribution, ...)
                       panel.qqmath(x, distribution = distribution, ...)
                   },
                   xlab = "log of Mean Expression",
                   ylab = "Quantiles of Standard Normal")
    }
    else if (model == "gamma")
    {
        ans <- 
            qqmath(~ means | mean.ranks, 
                   scales = list(relation = "free", draw = FALSE),
                   prepanel = function(x, ...) {
                       shape.gg <- (mean(x))^2 / var(x)
                       xx <- qgamma(ppoints(length(x)), shape = shape.gg)
                       list(xlim = range(xx), ylim = range(x))
                   },
                   panel = function(x, y, distribution, ...) {
                       shape.gg <- (mean(x))^2 / var(x)
                       xx <- qgamma(ppoints(length(x)), shape = shape.gg)
                       panel.qqmathline(x, distribution = function(p) qgamma(p, shape = shape.gg), ...)
                       panel.qqmath(x, distribution = function(p) qgamma(p, shape = shape.gg), ...)
                   },
                   xlab = "Mean Expression",
                   ylab = "Quantiles of fitted Gamma distribution")
    }
    ans
}














