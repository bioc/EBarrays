### Copyright 2003, Christina Kendziorski <kendzior@biostat.wisc.edu>,
### Michael Newton <newton@biostat.wisc.edu> and Deepayan Sarkar
### <deepayan@stat.wisc.edu>
###
### This file is part of the EBarrays library for R.  It is made
### available under the terms of the GNU General Public License,
### version 2, or at your option, any later version, incorporated
### herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


## slightly bad thing about this: bounds for optim are constrained






eb.createFamilyGG <-
    function()
{

    f0 <- function(theta, args)
    {
        ## alpha=theta[1], alpha0=theta[2], nu=theta[3]
        theta[2] * log(theta[3]) +
            lgamma(args$n * theta[1] + theta[2]) -
                args$n * lgamma(theta[1]) - lgamma(theta[2]) + 
                    ( (theta[1] - 1) * args$lprod.data) -
                        ( (args$n * theta[1] + theta[2]) * log(args$sum.data + theta[3])  )
    }


    f0.pp <- function(theta, args)
    {
        ## alpha=theta[1], alpha0=theta[2], nu=theta[3]
        ll1 <- theta[2] * log(theta[3]) +
            lgamma(args$n * theta[1] + theta[2]) -
                args$n * lgamma(theta[1]) - lgamma(theta[2])

        ll2 <- - (args$n * theta[1] + theta[2]) * (log(args$sum.data + theta[3]))

        ll1 + ll2
    }


    ## x is log intensity. Calculating marginal density of x
    logDensity <- function(x, theta)
    {

        ## returns log density of natural log of intensity
        ## theta = c(alpha, alpha0, nu)

        ## NB: separate terms for (alpha0 - 1) log(exp(x)) + x cancel out partially 

        lgamma(theta[1] + theta[2]) - lgamma(theta[1]) -    # scalar
            lgamma(theta[2]) + theta[2] * log(theta[3]) +   # scalar
                theta[1] * x - (theta[1] + theta[2]) * log(theta[3] + exp(x))
    }


    f0.arglist <- function(data, patterns)
    {

        ## returns a list with two components, common.args and
        ## pattern.args. common.args is a list of arguments to f0 that
        ## don't change from one pattern to another, whereas
        ## pattern.args[[i]][[j]] is a similar list of arguments, but
        ## specific to the columns in pattern[[i]][[j]]

        ## Note: f0 will also have an argument theta, which would need
        ## to be specified separately

        common.args <- list() ## nothing for GG
        pattern.args <- vector("list", length(patterns))
        for (i in seq(along = patterns))
        {
            pattern.args[[i]] <- vector("list", length(patterns[[i]]))
            for (j in seq(along = patterns[[i]]))
            {
                tmpdata <- data[, patterns[[i]][[j]], drop = FALSE]
                nn <- length(patterns[[i]][[j]])
                lprod.data <- rowSums(log(tmpdata))
                sum.data <- rowSums(tmpdata)
                pattern.args[[i]][[j]] <-
                    list(n = nn, lprod.data = lprod.data,
                         sum.data = sum.data)
            }
        }
        list(common.args = common.args, pattern.args = pattern.args)
    }
    


    thetaInit <- function(xx)
    {

        ## Function to take replicated array data xx and derive method of
        ## moment estimates for the shape parameters of the observation
        ## component and for the parameters of the inverse gamma mean
        ## component (leaves only the mixing proportions to be estimated)

        ## Based on the following theory:
        ## If X is an average expression value, according to the GG model,
        ##
        ## E[ X^q ] = numer/denom
        ##
        ##   where numer = nu^q Gamma( m*a+q) Gamma( a0-q )
        ##   and   denom = Gamma(m*a) Gamma(a0) m^q
        ## 
        ##  at least when q < a0, and where m is the number of replicates making X 



        ## cv <- function(x) { return(sqrt(var(x))/mean(x)) }
        cv <- function(x) { return(var(x)/(mean(x))^2) }

        tgamma <- function(x, quad)
        {
            tmp <- trigamma(x) - quad
            tmp
        } 

        ## ? xx <- log(as.matrix(xx))
        nrepx <- ncol(xx)

        ## now average within gene
        mx <- rowMeans(xx)

        ##    cv1 <- ( apply(xx, 1, cv) )^2     ## closer to unbiased to mean the squares
        cv1 <- ( apply(xx, 1, cv) )  ## closer to unbiased to mean the squares
        ax <- 1/(  mean(cv1) )
        ax <- 1/sqrt( mean(cv1) )

        ## so the above gets us estimates of the observation component shapes.

        qsupp <- seq(from = .001, to = 5, length=50 )
        bar.x <- outer(mx, qsupp, FUN="^")
        emp.x <- log( colMeans(bar.x) )

        qsupp2 <- qsupp^2

        fit.x <- lm( emp.x ~ qsupp + qsupp2 - 1 )  ##quadratic, no intercept

        quadCoef <- coef(fit.x)[2]
        if (quadCoef < .0005) stop("couldn't get reliable initial estimates")
        foo.x <- uniroot(tgamma, c(.001, 2500), quad = quadCoef)

        a0.x <- foo.x$root
        nu.x <- nrepx * exp( coef(fit.x)[1] + digamma(a0.x) - digamma(nrepx*ax) )
        ##x0.x <- ax * nu.x / a0.x

        c(alpha = ax, alpha0 = a0.x, nu = nu.x)

    }


    new("ebarraysFamily",
        "GG",
        description = "Gamma-Gamma",
        thetaInit = thetaInit,
        link = function(x) x, ## theta = (alpha, alpha_0, nu)
        invlink = function (x) {
            names(x) <- c("alpha", "alpha0", "nu")
            x
        },
        f0 = f0,
        f0.pp = f0.pp,
        logDensity = logDensity,
        f0.arglist = f0.arglist,
        lower.bound  = c(1.01, 0.01, 0.01),
        upper.bound = c(10000, 10000, 10000))
}







