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

checkModel <- function(data, fit, model = c("gamma", "lognormal", "lnnmv"),
                       number = 9, nb = 10, cluster=1, groupid=NULL){
  if (!require(lattice, quietly = TRUE)) stop("The lattice package could not be loaded")
  model <- match.arg(model)
  if (is(data, "exprSet") || is(data, "ExpressionSet")) data <- exprs(data)
  
  flag <- apply(data<=0,1,sum)
  data <- data[flag==0,]

  if(model=="gamma" | model=="lognormal"){
    clusters <- postprob(fit,data)$cluster
  
    pp.to.cluster <- function(pps){
      return(order(-pps)[1])
    }
  
    clusters <- apply(clusters,1,pp.to.cluster)
    data <- data[, unlist(fit@hypotheses@patterns[[1]])]
  
    curr.data <- data[clusters==cluster,]
  } else if (model=="lnnmv"){
    curr.data <- data[, unlist(fit@hypotheses@patterns[[1]])]
  }

  if (model == "gamma"){
    means <- apply(curr.data, 1, mean)
    mean.ranks <- equal.count(rank(means), number = number,
                              overlap = (number - length(means) / nb) / (number - 1))
    theta <- fit@family@invlink(fit@thetaEst[cluster,])
    ans <- qqmath(~ as.vector(curr.data) | mean.ranks, 
                  scales = list(relation = "free", draw = TRUE),
                  prepanel = function(x, ...) {
                    scale.gg <- mean(x) / theta[1]
                    xx <- qgamma(ppoints(length(x)), shape = theta[1], scale = scale.gg)
                    list(xlim = range(xx), ylim = range(x))
                  },
                  panel = function(x, ...) {
                    scale.gg <- mean(x) / theta[1]
                    xx <- qgamma(ppoints(length(x)), shape = theta[1], scale = scale.gg)
                    panel.qqmathline(x, distribution = function(p) qgamma(p, shape = theta[1], scale = scale.gg))
                    panel.qqmath(x, distribution = function(p) qgamma(p, shape = theta[1], scale = scale.gg))
                  },
                  xlab = "Quantiles of fitted Gamma distribution",
                  ylab = "Expression")
  } else if (model == "lognormal"){
    curr.data <- log(curr.data)
    means <- apply(curr.data, 1, mean)
    curr.data <- curr.data-outer(means, rep(1,dim(curr.data)[2]))
    mean.ranks <- equal.count(rank(means), number = number,
                              overlap = (number - length(means) / nb) / (number - 1))
    theta <- fit@family@invlink(fit@thetaEst[cluster,])
    ans <- qqmath(~ as.vector(curr.data) | mean.ranks,
                  distribution = function(p){qnorm(p, mean=0, sd=sqrt(theta[2]))},
                  scales = list(relation = "free", draw = TRUE),
                  panel = function(x, distribution, ...) {
                      panel.qqmathline(x, distribution = distribution, ...)
                      panel.qqmath(x, distribution = distribution, ...)
                    },
                  xlab = paste("Quantiles of Normal N(0,", round(theta[2],3),")", sep=""),
                  ylab = "log Expression centered at 0")
  } else if (model== "lnnmv"){
    curr.data <- log(curr.data)
    means <- apply(curr.data, 1, mean)
    if(is.null(groupid)) stop("Please specify groupid")
    vars <- fit@family@f0.arglist(data, fit@hypotheses@patterns, groupid)$common.args$sigmag2
    curr.data <- (curr.data-outer(means, rep(1,dim(curr.data)[2])))/outer(sqrt(vars), rep(1,dim(curr.data)[2]))
    mean.ranks <- equal.count(rank(means), number = number,
                              overlap = (number - length(means) / nb) / (number - 1))
    ans <- qqmath(~ as.vector(curr.data) | mean.ranks,
                  distribution = qnorm,
                  scales = list(relation = "free", draw = TRUE),
                  panel = function(x, distribution, ...) {
                      panel.qqmathline(x, distribution = distribution, ...)
                      panel.qqmath(x, distribution = distribution, ...)
                    },
                  xlab = "Quantiles of Standard Normal",
                  ylab = "Standardized log Expression")                
  }
  ans
}

checkVarsQQ <- function(data, groupid, ...){
  if (!require(lattice, quietly = TRUE)) stop("The lattice package could not be loaded")
  if (is(data, "exprSet") || is(data, "ExpressionSet")) data <- exprs(data)
  flag <- apply(data<=0,1,sum)
  data <- data[flag==0,]
  data <- log(data)
  n <- sum(groupid!=0)
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

  # gene specific sample variance
  svars <- vars/(n-ngroup)
  
  # The following is the moment estimate of the (v0,sigma02)
  v0 <- mean(svars)^2/((length(svars)-1)/length(svars)*var(svars))*2+4
  sigma02 <- mean(svars)*(v0-2)/v0
                  
  ans <- qqmath(~ svars, 
                distribution = function(p){v0*sigma02/qchisq(1-p, v0)},
                panel = function(x, distribution, ...) {
                  panel.qqmathline(x, distribution = distribution, ...)
                  panel.qqmath(x, distribution = distribution, ...)
                },               
                xlab = expression(paste("Quantiles of Inv-",
                  chi^2, "(", list(hat(nu)[0], hat(sigma)[0]^2), ")")),
                ylab = "Gene Specific Sample Variances", ...)
  ans
}

checkVarsMar <- function (data, groupid, xlab, ylab, ...){
  if (!require(lattice, quietly = TRUE)) stop("The lattice package could not be loaded")
  if (is(data, "exprSet") || is(data, "ExpressionSet")) data <- exprs(data)

  # density of scaled inverse chi-square distribution
  dsichisq <- function(x, df, scale){
    return(dchisq(df*scale/x, df)*df*scale/(x^2))
  }

  flag <- apply(data<=0,1,sum)
  data <- data[flag==0,]
  data <- log(data)
  nn <- sum(groupid!=0)
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
  svars <- vars/(nn-ngroup)

  # The following is the moment estimate of the (v0,sigma02)
  v0 <- mean(svars)^2/((length(svars)-1)/length(svars)*var(svars))*2+4
  sigma02 <- mean(svars)*(v0-2)/v0
  
  if (missing(xlab)) xlab <- "Gene Specific Sample Variances"
  if (missing(ylab)) ylab <- "Density"
 
  ans <- histogram( ~ svars, endpoints=c(min(svars), max(svars)),
                   type = "density", ...)
  y.limit <- ans$y.limit
  ans <- histogram( ~ svars, endpoints=c(min(svars), max(svars)),
                   xlab = xlab, ylab = ylab, type = "density",
                   prepanel = function(x,...){
                     xx <- seq(min(x), max(x), 0.0001)
                     yy <- dsichisq(xx,v0,sigma02)
                     list(ylim = c(range(yy)[1], max(range(yy)[2], y.limit[2])))
                   },
                   panel = function(x, ...) {
                     panel.histogram(x, ...)
                     xx <- seq(min(x), max(x), 0.0001)
                     yy <- dsichisq(xx,v0,sigma02)
                     panel.xyplot(xx, yy, type="l", col="red", lwd=2)
                   }, ... )
  ans
}

#plots
plotMarginal <- function (fit, data, kernel = "rect", n = 100, bw = "nrd0", adjust = 1, xlab, ylab, ...)
{

    ## Input: data is a matrix or exprSet. fit is the
    ## output from emfit

    ## Output: Plot comparing the marginal distribution under the
    ## theoretical model to the empirical distribution, on the log
    ## scale.

    ## Note: x_ij values have the same marginal, but x_i1, ... x_im
    ## are not necessarily independent. Still OK to estimate density
    ## from all data, since the result can be thought of as the
    ## average of several density estimates, one for each column.

    if (!require(lattice, quietly = TRUE)) stop("The lattice package could not be loaded")

    n.cluster<-dim(fit@thetaEst)[1]

    if (is(data, "exprSet") || is(data, "ExpressionSet")) data <- exprs(data)

    flag <- apply(data<=0,1,sum)
    data <- data[flag==0,]

    clusters <- postprob(fit,data)$cluster

    pp.to.cluster<-function(pps)
    {
        return(order(-pps)[1])
    }

    clusters<-apply(clusters,1,pp.to.cluster)


    data <- log(data)

    #lims <- range(data)

    ## By default, plot average shifted histogram of log expressions
    ## instead of plain histogram. Other kernels also possible

    #logmarg<-NULL
    #y<-NULL
    #supp<-NULL
    #cluster<-NULL
    data <- data[,unlist(fit@hypotheses@patterns[[1]])]
    empmarg <- density(data, bw = "nrd0", adjust = 1,
                       kernel = kernel, n = n)

    all.marg<-0
    ps<-apply(fit@probEst,1,sum)
    
    for (i in 1:n.cluster)
    {
        theta <- fit@family@invlink(fit@thetaEst[i,])
        all.marg <- all.marg+ps[i]*exp(fit@family@logDensity(empmarg$x, theta))
    }
        
    supp <- empmarg$x
    marg <- all.marg
    y <- empmarg$y
    cluster <- rep(paste("All Genes:",dim(data)[1],"genes"),length(empmarg$x))

    if (n.cluster>1) {
        for (i in 1:n.cluster) {
            curr.data<-data[clusters==i,]
        
            if (sum(curr.data>0)>100)
            {
                empmarg <- density(curr.data, bw = "nrd0", adjust = 1,
                                   kernel = kernel, n = n)
                                        #from = lims[1], to = lims[2])
                
                supp <- c(supp, empmarg$x)
                theta <- fit@family@invlink(fit@thetaEst[i,])
                marg <- c(marg,exp(fit@family@logDensity(empmarg$x, theta)))
                
                if (i>9)
                {
                    cluster<-c(cluster,rep(paste("Cluster",i,":",sum(clusters==i),"genes"),length(empmarg$x)))
                } else {
                    cluster<-c(cluster,rep(paste("Cluster ",i,":",sum(clusters==i),"genes"),length(empmarg$x)))
                }
            
                y<-c(y,empmarg$y)
            }
        }
    }

    if (missing(xlab)) xlab <- "log(expressions)"
    if (missing(ylab)) ylab <- "Density"

    xyplot(marg + y ~ supp|cluster, type = 'l', as.table=TRUE,
           allow.mult = TRUE, scales=list(y="free"),
           key = simpleKey(text = c("Marginal Density of log Expressions from Fitted Model",
                           "Empirical Kernel Density of log Expressions"),
           lines = TRUE, points = FALSE),
           xlab = xlab,
           ylab = ylab, ...)
}


plotCluster<-function (fit, data, cond = NULL, ncolors = 123, sep=TRUE, transform=NULL) 
{
    pps<-postprob(fit,data)
    clu<-pps$cluster
    pat<-pps$pattern[,2]

    pp.to.clus<-function(pps)
    {
        return(order(-pps)[1])
    }

    clus<-apply(clu,1,pp.to.clus)
    
    o<-order(clus*10+pat)
    z<-(log(data[o,]))
    
    if (!is.null(cond)) {
      o<-order(cond)
      z<-t(z[,o])
    }
    
    zr <- range(z, na.rm = TRUE)
    zmax <- max(abs(zr))
    zmin <- zr[1]

    low <- c(0, 1, 0)
    high <- c(1, 0, 0)

    zerocenter <- (zmin < 0)
    if (zerocenter) {
        zlim <- c(-zmax, zmax)
    } else {
        zlim <- c(zmin, zmax)
    }

    col <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
        high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
    x<-1:nrow(z)
    y<-1:ncol(z)

    if (!is.null(transform)) {
      image(y,x,t(transform(z)), zlim = transform(zlim), col = col, axes = FALSE, xlab = "Replicate", ylab = "Gene")
    } else {
      image(y,x,t(z), zlim = zlim, col = col, axes = FALSE, xlab = "Replicate", ylab = "Gene")
    }
    #axis(1, at = x)
    #len<-max(10,length(x)+5)
    #axis(2, at = seq(1,floor(ncol(z)/len)*len,length=len))

    if (sep==TRUE) {
      clus<-table(clus)
    
      clusline<-0
      for (i in 1:(length(clus)-1))
      {
        clusline<-clusline+clus[i]
        abline(h=clusline,lwd=3)
      }
    
      if (!is.null(cond)) {
        cond<-table(cond)
    
        condline<-0.5
        for (i in 1:(length(cond)-1))
          {
            condline<-condline+cond[i]
            abline(v=condline)
          }
      }
    }
}

plot.ebarraysEMfit<-function(x,data,plottype="cluster",...)
{
    if (plottype=="cluster")
    {
        plotCluster(fit=x,data=data,...)
    } else if (plottype=="marginal") {
        plotMarginal(fit=x,data=data,...)
    }
}
