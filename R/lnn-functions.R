
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


## Derivation:

## Given matrix X (kxn) with k genes and n rows, where x_1, ..., x_k
## are rows of X. Each x_i is iid

##        MVN(\mu 1_n, \sigma^2 I_n + \tao^2 M_n)

## where M_n = 11' = nxn matrix of ones. Then, with

## A1 = 1/k \Sum_i [ x_i^T x_i ]
## A2 = 1/k \Sum_i [ 1' x_i ]^2
## A3 = 1/k \Sum_i [ 1' x_i ]

## the log likelihood reduces to

##  l(sigma^2, tao^2, mu) = - n/2 log(2 pi) - (n-1)/2 log(sigma^2)
##  - 1/2 log(sigma^2 + n tao^2) - c/2 (n mu^2 + A1 - 2 mu A3)
##  - d/2 (A2 + n^2 mu^2 - 2n mu A3) 

## where c = 1/sigma^2, d = - tao^2 / (sigma^2 (sigma^2 + n tao^2))







eb.createFamilyLNN <-
    function()
{

    f0 <- function(theta, args)
    {
        ## mu0=theta[1], log(sigma2)=theta[2], log(tau02)=theta[3]

        theta[2:3] <- exp(theta[2:3])
        spnt <- theta[2] + args$n * theta[3] ## sigma^2 + n tao0^2
        negloglik <-

            c(args$n, args$n - 1, 1) %*%
                log(c(2 * base::pi, theta[2], spnt)) +
                    
                    (args$n * theta[1] * theta[1] + args$sumsq.logdata - 2
                     * theta[1] * args$sum.logdata) / theta[2] -
                         
                         theta[3] * (args$sum.logdata * args$sum.logdata +
                                     args$n * args$n * theta[1] * theta[1]
                                     - 2 * args$n * theta[1] *
                                     args$sum.logdata ) / (theta[2] *
                                                           spnt)

        ##print(str(negloglik))
        -0.5 * negloglik
    }


    ## x is log intensity. Calculating marginal density of x
    logDensity <- function(x, theta)
    {
        ## returns log density of natural log of intensity
        ## theta = c(mu0, sigma^2, tao0^2)

        dnorm(x, mean = theta[1],
              sd = sqrt(sum(theta[2:3])),
              log = TRUE)
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

        data <- log(data)
        k <- nrow(data)
        common.args <- list() ## nothing for LNN
        pattern.args <- vector("list", length(patterns))
        for (i in seq(along = patterns))
        {
            pattern.args[[i]] <- vector("list", length(patterns[[i]]))
            for (j in seq(along = patterns[[i]]))
            {
                tmpdata <- data[, patterns[[i]][[j]], drop = FALSE]
                pattern.args[[i]][[j]] <-
                    list(n = length(patterns[[i]][[j]]),
                         sum.logdata = rowSums(tmpdata),
                         sumsq.logdata = rowSums(tmpdata * tmpdata))
            }
        }
        list(common.args = common.args, pattern.args = pattern.args)
    }


    ## function to get initial values of theta if they are not supplied
    thetaInit <- function(data) {
        data <- log(as.matrix(data))
        ncols <- ncol(data)
        means <- rowMeans(data) 
        vars <- apply(data, 1, FUN=var)
        sigmasq.hat <- mean(vars)
        mu.hat <- mean(means)
        tausq.hat <- var(means)
        ## return value: theta = (mu_0, log(sigma^2), log(tao_0^2))
        c(mu.hat, log(sigmasq.hat), log(tausq.hat))
    }


    
    new("ebarraysFamily",
        "LNN",
        description = "Lognormal-Normal",
        thetaInit = thetaInit,
        link = function(x) {
            x[2:3] <- log(x[2:3])
            ## theta = (mu_0, log(sigma^2), log(tao_0^2))
            x 
        },
        invlink = function (x) {
            x[2:3] <- exp(x[2:3])
            names(x) <- c("mu_0", "sigma^2", "tao_0^2")
            x
        },
        f0 = f0,
        f0.pp = f0,
        logDensity = logDensity,
        f0.arglist = f0.arglist,
        lower.bound = c(-Inf, -Inf, -Inf), 
        upper.bound = c(Inf, Inf, Inf))
}




# eb.createFamilyLNN.old <-
#     function()
# {

#     f0 <- function(theta, args)
#     {
#         ## mu0=theta[1], sigma2=theta[2], tau02=theta[3]

#         spnt <- theta[2] + args$n * theta[3] ## sigma^2 + n tao0^2
#         negloglik <-

#             c(args$n, args$n - 1, 1) %*%
#                 log(c(2 * base::pi, theta[2], spnt)) +
                    
#                     (args$n * theta[1] * theta[1] + args$sumsq.logdata - 2
#                      * theta[1] * args$sum.logdata) / theta[2] -
                         
#                          theta[3] * (args$sum.logdata * args$sum.logdata +
#                                      args$n * args$n * theta[1] * theta[1]
#                                      - 2 * args$n * theta[1] *
#                                      args$sum.logdata ) / (theta[2] *
#                                                            spnt)

#         ##print(str(negloglik))
#         -negloglik
#     }

#     f0.arglist <- function(data, patterns)
#     {

#         ## returns a list with two components, common.args and
#         ## pattern.args. common.args is a list of arguments to f0 that
#         ## don't change from one pattern to another, whereas
#         ## pattern.args[[i]][[j]] is a similar list of arguments, but
#         ## specific to the columns in pattern[[i]][[j]]

#         ## Note: f0 will also have an argument theta, which would need
#         ## to be specified separately

#         data <- log(data)
#         k <- nrow(data)
#         common.args <- list() ## nothing for LNN
#         pattern.args <- vector("list", length(patterns))
#         for (i in seq(along = patterns))
#         {
#             pattern.args[[i]] <- vector("list", length(patterns[[i]]))
#             for (j in seq(along = patterns[[i]]))
#             {
#                 tmpdata <- data[, patterns[[i]][[j]], drop = FALSE]
#                 pattern.args[[i]][[j]] <-
#                     list(n = length(patterns[[i]][[j]]),
#                          sum.logdata = rowSums(tmpdata),
#                          sumsq.logdata = rowSums(tmpdata * tmpdata))
#             }
#         }
#         list(common.args = common.args, pattern.args = pattern.args)
#     }

#     new("ebarraysFamily",
#         "LNN",
#         description = "Lognormal-Normal",
#         link = function(x) x, ## theta = (mu_0, sigma^2, tao_0^2)
#         invlink = function (x) {
#             names(x) <- c("mu_0", "sigma^2", "tao_0^2")
#             x
#         },
#         f0 = f0,
#         f0.pp = f0,
#         f0.arglist = f0.arglist,
#         lower.bound  = c(-10000, 0.01, 0.01),
#         upper.bound = c(10000, 10000, 10000))
# }







# f0lnn.pp <-
#     function(theta, args)
# {
#     ## mu0=theta[1], sigma2=theta[2], tau02=theta[3]

#     spnt <- theta[2] + args$n * theta[3] ## sigma^2 + n tao0^2
#     negloglik <-

#         c(args$n, args$n - 1, 1) %*%
#             log(c(2 * base::pi, theta[2], spnt)) +

#                 (args$n * theta[1] * theta[1] + args$sumsq.logdata - 2
#                  * theta[1] * args$sum.logdata) / theta[2] -

#                      theta[3] * (args$sum.logdata * args$sum.logdata +
#                                  args$n * args$n * theta[1] * theta[1]
#                                  - 2 * args$n * theta[1] *
#                                  args$sum.logdata ) / (theta[2] *
#                                                        spnt)

#     -negloglik
# }

