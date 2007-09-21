setClass("ebarraysFamily",
         representation("character",
                        description = "character",
                        thetaInit = "function",
                        link = "function",  ## theta = f(model par)
                        invlink = "function", ## inverse
                        f0 = "function",
                        f0.pp = "function",
                        logDensity = "function",
                        f0.arglist = "function",
                        lower.bound = "numeric",
                        upper.bound = "numeric",
                        deflate = "function"))

setAs("character", "ebarraysFamily",

      function(from) {

          if (from == "GG") 
              eb.createFamilyGG()
          else if (from == "LNN")
              eb.createFamilyLNN()
          else if (from == "LNNMV")
              eb.createFamilyLNNMV()
          else stop(paste("\n The only families recognized by name are GG (Gamma-Gamma)",
                          "\n and LNN (Lognormal-Normal). Other families need to be",
                          "\n specified as objects of class 'ebarraysFamily'"))
          
      })

setMethod("show", "ebarraysFamily",

          function(object) {
              cat(paste("\n Family:", object, "(", object@description, ")",
                        "\n\nAdditional slots:\n"))
              show(getSlots(class(object)))
              invisible(object)
          })

##ebarraysFamily

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

    deflate <- function(theta, args, num.sample=10000)
    {
        ## alpha=theta[1], alpha0=theta[2], nu=theta[3]

        num.groups<-length(args)
        a <- rep(0,num.groups)
        b <- matrix(0,length(args[[1]]$sum.data),num.groups)
    
        for (i in 1:num.groups)
        {
          a[i]<-args[[i]]$n*theta[1]+theta[2]
          b[,i]<-args[[i]]$sum.data+theta[3]
        }
        
        if (num.groups==2)
        {
          ans <- pbeta(b[,1]/(b[,1]+b[,2]),a[1],a[2])
        } else {
          V.sample<-matrix(0,num.sample,num.groups)
          for (i in 1:num.groups)
          {
            V.sample[,i]<-rgamma(num.sample,shape=a[i],rate=1)
          }

          ans <- rep(0,length(args[[i]]$sum.data))
          
          for (i in 1:length(args[[i]]$sum.data))
          {
            ans.vec<-rep(TRUE,num.sample)
            for (j in 1:(num.groups-1))
            {
              ans.vec<-ans.vec & (V.sample[,j]/b[i,j]>V.sample[,j+1]/b[i,j+1])
            }
            ans[i]<-mean(ans.vec)
          }
        }
        return(log(gamma(num.groups+1)*ans))
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
        upper.bound = c(10000, 10000, 10000),
        deflate = deflate)
}

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

    deflate <- function(theta, args, num.sample=10000)
    {
        ## mu0=theta[1], log(sigma2)=theta[2], log(tau02)=theta[3]

        theta[2:3] <- exp(theta[2:3])
        num.groups<-length(args)
        spnt <- rep(0,num.groups)
        a <- matrix(0,length(args[[1]]$sum.logdata),num.groups)
        b <- rep(0,num.groups)
    
        for (i in 1:num.groups)
        {
          spnt[i] <- theta[2] + args[[i]]$n * theta[3] ## sigma^2 + n tao0^2        
          a[,i]<-(theta[2]*theta[1] + theta[3]*args[[i]]$sum.logdata)/spnt[i]
          b[i]<-theta[2] * theta[3]/spnt[i]
        }
        
        if (num.groups==2)
        {
          ans <- pnorm((a[,1]-a[,2])/sqrt(b[1]^2+b[2]^2))
        } else {
          V.sample<-matrix(rnorm(num.sample*num.groups),num.sample,num.groups)
          for (i in 1:num.groups)
          {
            V.sample[,i]<-V.sample[,i]*b[i]
          }

          ans <- rep(0,length(args[[i]]$sum.logdata))
          
          for (i in 1:length(args[[i]]$sum.logdata))
          {
            ans.vec<-rep(TRUE,num.sample)
            for (j in 1:(num.groups-1))
            {
              ans.vec<-ans.vec & (V.sample[,j]+a[i,j]>V.sample[,j+1]+a[i,j+1])
            }
            ans[i]<-mean(ans.vec)
          }
        }
 
        return(log(gamma(num.groups+1)*ans))
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
        upper.bound = c(Inf, Inf, Inf),
        deflate = deflate)
}
