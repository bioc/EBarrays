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

##        MVN(\mu 1_n, \sigma_g^2 I_n + \tao^2 M_n)

## where M_n = 11' = nxn matrix of ones. Then, with

## A1 = x_i x_i^T
## A2 = [ x_i 1' ]^2
## A3 = [ x_i 1' ]

## the log likelihood reduces to

##  l(sigma_g^2, tao^2, mu) = - n/2 log(2 pi) - (n-1)/2 log(sigma_g^2)
##  - 1/2 log(sigma_g^2 + n tao^2) - c/2 (n mu^2 + A1 - 2 mu A3)
##  - d/2 (A2 + n^2 mu^2 - 2n mu A3) 

## where c = 1/sigma_g^2, d = - tao^2 / (sigma_g^2 (sigma_g^2 + n tao^2))

eb.createFamilyLNNMV <-
    function()
{

    f0 <- function(theta, args)
    {
        ## mu0=theta[1], log(tau02)=theta[2]

        theta[2] <- exp(theta[2])
        spnt <- args$sigmag2 + args$n * theta[2] ## sigma_g^2 + n tao0^2
        negloglik <-

            log(cbind(rep(2 * base::pi,length(args$sigmag2)), args$sigmag2, spnt)) %*%
                c(args$n, args$n - 1, 1)+
                    
                    (args$n * theta[1] * theta[1] + args$sumsq.logdata - 2
                     * theta[1] * args$sum.logdata) / args$sigmag2 -
                         
                         theta[2] * (args$sum.logdata * args$sum.logdata +
                                     args$n * args$n * theta[1] * theta[1]
                                     - 2 * args$n * theta[1] *
                                     args$sum.logdata ) / (args$sigmag2 *
                                                           spnt)

        ##print(str(negloglik))
        -0.5 * negloglik
    }

    
    ## x is log intensity. Calculating marginal density of x
    logDensity <- function()
    {
      cat("\n logDensity is not available for family \"LNNMV\" \n")
    }


    f0.arglist <- function(data, patterns, groupid)
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
        n <- sum(groupid!=0)
        
        ## Assume sigma_g^2 ~ Inv-chi-square (v0,sigma0^2).
        groups <- unique(groupid)
        groups <- groups[groups!=0]
        ngroup <- length(groups)
        vars <- 0
        for(i in 1:ngroup){
          temp <- data[,groupid==groups[i]]
          ncoltemp <- sum(groupid==groups[i])
          if(ncoltemp==1){
              vars <- vars  
          }else{
              vars <- vars + (ncoltemp-1)*apply(temp,1,FUN=var)
          }
        }
        ## gene specific sample variance
        svars <- vars/(n-ngroup)

        # The following is the moment estimate of the (v0,sigma02)
        v0 <- mean(svars)^2/((length(svars)-1)/length(svars)*var(svars))*2+4
        sigma02 <- mean(svars)*(v0-2)/v0
                
        ## estimate of sigma_g^2 for each gene, do not depend on patterns
        common.args <- list(sigmag2=(v0*sigma02+vars)/(n+v0-2))
        
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
        means <- rowMeans(data) 
        mu.hat <- mean(means)
        tausq.hat <- var(means)
        ## return value: theta = (mu_0, log(tao_0^2))
        c(mu.hat, log(tausq.hat))
    }

    deflate <- function(theta, args, num.sample=10000)
    {
        stop("LNNMV family does not support ordered hypotheses.")
    }
    
    new("ebarraysFamily",
        "LNNMV",
        description = "Lognormal-Normal with modified variances",
        thetaInit = thetaInit,
        link = function(x) {
            x[2] <- log(x[2])
            ## theta = (mu_0, log(tao_0^2))
            x 
        },
        invlink = function (x) {
            x[2] <- exp(x[2])
            names(x) <- c("mu_0", "tao_0^2")
            x
        },
        f0 = f0,
        deflate = deflate,
        f0.pp = f0,
        logDensity = logDensity,
        f0.arglist = f0.arglist,
        lower.bound = c(-Inf, -Inf), 
        upper.bound = c(Inf, Inf))
}
