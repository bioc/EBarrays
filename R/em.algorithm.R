###class for result of single time emfit
setClass("ebarraysEMfit",
         representation(family = "ebarraysFamily",
                        hypotheses = "ebarraysPatterns",
                        thetaEst = "matrix",
                        probEst = "matrix",
                        ic="numeric"))

setMethod("show", "ebarraysEMfit",
          function(object) {

              cat(paste("\n EB model fit",
                        "\n\t Family:", object@family,
                        "(", object@family@description, ")"))
              n.cluster = dim(object@probEst)[1]
              cat("\n\n Model parameter estimates:\n\n")
              thetaEst<-NULL
              for (i in 1:dim(object@thetaEst)[1])
                  thetaEst<-rbind(thetaEst,object@family@invlink(object@thetaEst[i,]))

              if (is.null(rownames(object@thetaEst))) {
                rownames(thetaEst)<-paste("Cluster",1:dim(thetaEst)[1])
              }
              
              print(data.frame(thetaEst))
              cat("\n Estimated mixing proportions:\n\n")
              probEst<-object@probEst
              colnames(probEst)<-paste("Pattern",1:dim(probEst)[2])
              if (is.null(rownames(object@probEst))) {
                rownames(probEst)<-paste("Cluster",1:dim(thetaEst)[1])
              }
              print(data.frame(probEst))
          })

setGeneric("emfit",
           function(data, family, hypotheses, ...) standardGeneric("emfit"))


setMethod("emfit",
          signature(data = "matrix",
                    family = "character",
                    hypotheses = "ebarraysPatterns"),

          function(data,
                   family,
                   hypotheses, ...) {

              family <- as(family, "ebarraysFamily")
              callGeneric(data = data,
                          family = family,
                          hypotheses = hypotheses, ...)

          })



setMethod("emfit",
          signature(data = "matrix",
                    family = "ebarraysFamily",
                    hypotheses = "ebarraysPatterns"),

          function(data,
                   family,
                   hypotheses,
                   cluster=1,
                   type=2,
                   criterion="BIC",
                   cluster.init = NULL,
                   independent = F,
                   num.iter = 20,
                   verbose = getOption("verbose"),
                   optim.control = list(),
                   ...) {

             if (!require(cluster, quietly = TRUE)) stop("The cluster package could not be loaded")

             if (family=="LNNMV") {
                return(emfitmv(data=data,
                              family=family,
                              hypotheses=hypotheses,
                              num.iter=num.iter,
                              verbose=verbose,
                              optim.control=optim.control,
                              ...))
             } else if (type==1) { #fixed cluster
                return(emfit2(data,
                              family,
                              hypotheses,
                              cluster,
                              num.iter,
                              verbose,
                              optim.control,
                              ...))
              } else if (type==2) { # one number of clusters
                if (length(cluster)==1) {
                  return(emfit1(data,
                                family,
                                hypotheses,
                                cluster,
                                cluster.init,
                                independent,
                                num.iter,
                                verbose,
                                optim.control,
                                ...))
                } else { # several numbers of clusters
                  return(emfit3(data,
                                family,
                                hypotheses,
                                cluster,
                                criterion,
                                cluster.init,
                                independent,
                                num.iter,
                                verbose,
                                optim.control,
                                ...))
                }
              }
           })

emfit0 <- function(data,
                   family,
                   hypotheses,
                   num.iter = 20,
                   verbose = getOption("verbose"),
                   optim.control = list(),
                   ...) {

              family <- as(family, "ebarraysFamily")

              ## Store current process time
              if (exists("proc.time")) startTime <- proc.time()

              ## Define function to be optimized by optim:

              complete.loglik <-
                  function(theta,
                           zhat,
                           zhatlogp,
                           hypotheses,
                           f0, f0.args)
                  {
                      pred.dens <- zhat # values unimportant, just need the same dimensions

                      ## NOTE: f0 returns the log value of the predictive probability
                      ## density function.
                      for (i in seq(along = hypotheses))
                      {
                          foo <- 0 #numeric(nrow(zhat))
                          for (j in seq(along = hypotheses[[i]]))
                          {
                              foo <- 
                                  foo +
                                      f0(theta = theta,
                                         args = c(f0.args$common.args,
                                         f0.args$pattern.args[[i]][[j]]))
                          }
                          pred.dens[, i] <- foo
                      }
                      -(sum(zhat * pred.dens) + zhatlogp)
                  }


              ## Some preprocessing:
              verbose <- as.numeric(verbose)
              patternObj <- hypotheses
              hypotheses <- as(hypotheses@patterns, "list")

              if (verbose > 0) cat(paste("\n Checking for negative entries..."))
              all.pos <- apply(data, 1, function(x) !any(x <= 0))
              if (verbose > 0 && sum(all.pos) != length(all.pos)) {
                  cat(paste("\n\t", length(all.pos) - sum(all.pos),
                            "rows out of", length(all.pos),
                            "had at least one negative entry"))
                  cat("\n\t These rows will not be used in the EM fit\n")
              }


              if (verbose > 0)
                  cat(paste("\n Generating summary statistics for patterns.",
                            "\n This may take a few seconds...\n"))


              f0.args <-
                          family@f0.arglist(data[all.pos, ],
                                    hypotheses)


              theta.temp <- family@thetaInit(data[all.pos, ])
              ## should we rm(data) now ? Is data copied in memory ?


              p.init <- rep(1, length(hypotheses))
              p.init <- p.init / sum(p.init)
              p.temp <- p.init

              ## Start Iterations
              if (verbose > 0)
                  cat(paste("\n Starting EM iterations (total", num.iter,
                            ").\n This may take a while\n\n"))

              notdone <- TRUE ## not used now, perhaps some convergence criteria
              iter <- 0
              zhat <- matrix(0, sum(all.pos), length(hypotheses))

              while ({iter <- iter + 1;  (iter <= num.iter && notdone)})
              {

                  if (verbose > 0)
                      cat(paste("\t Starting iteration", iter, "...\n"))
                  
                  ## E step
                  zhat[,] <- 0

                  ## NOTE: f0 returns the log value of the predictive
                  ## probability density function. Thus to get the log of the
                  ## product of the predictive probability densities, we add
                  ## (log(PA=1 * PA=2) = log(PA=1) + log(PA=2)).

                  for (i in seq(along = hypotheses))
                  {
                      foo <- 0
                      for (j in seq(along = hypotheses[[i]]))
                      {
                          foo <- 
                              foo + 
                                  family@f0(theta = theta.temp,
                                            args = c(f0.args$common.args,
                                            f0.args$pattern.args[[i]][[j]]))
                      }
                      #zhat[, i] <- foo

                      if (patternObj@ordered & length(hypotheses[[i]])>1)
                      {
                        zhat[,i] <- foo+family@deflate(theta = theta.temp,f0.args$pattern.args[[i]])
                      } else {
                        zhat[,i] <- foo
                      }
                  }

                  zhat <- .Call("makeProbVect", zhat, p.temp, PACKAGE = "EBarrays")

                  #complete.lik<-zhat #value unimportant, just need the same dimension
                  #for (i in 1:dim(zhat)[1])
                  #{
                  #  complete.lik[i,]<-exp(zhat[i,])*p.temp
                  #}                  
                  #zhat<-complete.lik/apply(complete.lik,1,sum)

                  p.temp <- colSums(zhat)
                  p.temp <- p.temp / sum(p.temp)
                  if (verbose > 1) {
                      cat("\n\n      p.temp: ")
                  }


                  ## M step
                  zhatlogp <- sum(zhat %*% log(p.temp))
                  fit.temp <-
                      optim(par = theta.temp,
                            fn = complete.loglik,
                            method = c("L-BFGS-B"),
                            lower = family@lower.bound,
                            upper = family@upper.bound,
                            control = optim.control,
                            zhat = zhat,
                            zhatlogp = sum(zhat %*% log(p.temp)),
                            hypotheses = hypotheses,
                            f0 = family@f0,
                            f0.args = f0.args)
                  theta.temp <- fit.temp$par

                  if (verbose > 1) {
                      cat("\n\n  theta.temp: ")
                      print(theta.temp)
                  }
              }


              ## Finish up

              if (verbose > 1) {
                  cat("\n\n")
                  print("model parameter estimates:")
                  print(theta.temp)
                  print("Estimates of mixing proportions:")
                  print(p.temp)
              }

              if (verbose && exists("proc.time")) {
                  endTime <- proc.time()
                  cat(sprintf("\n\n Fit used %5.2f seconds user time\n",
                              endTime[1] - startTime[1]))
              }

              thetaEst<-rbind(NULL,theta.temp)
              probEst<-rbind(NULL,p.temp)
              rownames(thetaEst)<-NULL
              rownames(probEst)<-NULL

              new("ebarraysEMfit", 
                   family = family,
                   hypotheses = patternObj,
                   thetaEst = thetaEst,
                   probEst = probEst,
                   ic = fit.temp$value+(3+length(hypotheses)-1)*c(0.5*log(sum(all.pos)),1,log(log(sum(all.pos)))))
            }

emfit1<-function(data,
                 family,
                 hypotheses,
                 n.cluster,
                 cluster.init = NULL,
                 independent = FALSE,
                 num.iter = 20,
                 verbose = getOption("verbose"),
                 optim.control = list(),
                 ...) {

              if (is(data, "ExpressionSet")) data<-exprs(data)
              
              family <- as(family, "ebarraysFamily")
              
              if (!is.null(cluster.init) & length(cluster.init)!=dim(data)[1]) {
                warning("Number of genes appeared in cluster.init does not match the data; cluster.init is discarded!")
                cluster.init <- NULL
              }

              ## Store current process time
              if (exists("proc.time")) startTime <- proc.time()

              ## Define function to be optimized by optim:

              complete.loglik <-
                  function(theta,
                           zhat,
                           zhatlogp,
                           hypotheses,
                           f0, f0.args)
                  {
                      pred.dens <- zhat # values unimportant, just need the same dimensions

                      ## NOTE: f0 returns the log value of the predictive probability
                      ## density function.
                      for (i in seq(along = hypotheses))
                      {
                          foo <- 0 #numeric(nrow(zhat))
                          for (j in seq(along = hypotheses[[i]]))
                          {
                              foo <- 
                                  foo +
                                      f0(theta = theta,
                                         args = c(f0.args$common.args,
                                         f0.args$pattern.args[[i]][[j]]))
                          }
                          pred.dens[, i] <- foo
                      }
                      -(sum(zhat * pred.dens) + zhatlogp)
                  }


              ## Some preprocessing:
              verbose <- as.numeric(verbose)
              patternObj <- hypotheses
              hypotheses <- as(hypotheses@patterns, "list")

              if (verbose > 0) cat(paste("\n Checking for negative entries..."))
              all.pos <- apply(data, 1, function(x) !any(x <= 0))
              if (verbose > 0 && sum(all.pos) != length(all.pos)) {
                  cat(paste("\n\t", length(all.pos) - sum(all.pos),
                            "rows out of", length(all.pos),
                            "had at least one negative entry"))
                  cat("\n\t These rows will not be used in the EM fit\n")
              }

              data<-data[all.pos,]

              if (is.null(cluster.init)) {
                  cluster.init<-clara(log(data),n.cluster)$clustering
              } else {
                  cluster.init<-cluster.init[all.pos,]
              }

              unique.cluster<-unique(cluster.init)
              n.cluster <- length(unique.cluster)
              n.gene<-dim(data)[1]
              
              if (n.cluster==1) {
                return(emfit0(data,
                              family,
                              patternObj,
                              num.iter,
                              verbose,
                              optim.control,
                              ...))
              }
                   
              if (verbose > 0)
                  cat(paste("\n Generating summary statistics for patterns.",
                            "\n This may take a few seconds...\n"))

              f0.args <- family@f0.arglist(data,hypotheses)

              theta.temp<-NULL
              for (i in unique(cluster.init))
              {
                  theta.temp<-rbind(theta.temp,family@thetaInit(data[cluster.init==i,]))
              }
              
              p.temp <- matrix(1, n.cluster,length(hypotheses))

              p.temp <- p.temp / sum(p.temp)

              ## Start Iterations
              if (verbose > 0)
                  cat(paste("\n Starting EM iterations (total", num.iter,
                            ").\n This may take a while\n\n"))

              notdone <- TRUE ## not used now, perhaps some convergence criteria
              iter <- 0
              zhat<-array(0,dim=c(n.gene, n.cluster, length(hypotheses)))
              
              while ({iter <- iter + 1;  (iter <= num.iter && notdone)})
              {

                  if (verbose > 0)
                      cat(paste("\t Starting iteration", iter, "...\n"))
                  
                  ## E step
                  zhat<-zhat*0
                  
                  ## NOTE: f0 returns the log value of the predictive
                  ## probability density function. Thus to get the log of the
                  ## product of the predictive probability densities, we add
                  ## (log(PA=1 * PA=2) = log(PA=1) + log(PA=2)).

                  for (k in 1:n.cluster)
                  {
                    for (i in seq(along = hypotheses))
                    {
                      foo <- 0
                      for (j in seq(along = hypotheses[[i]]))
                      {
                          foo <- 
                              foo + 
                                  family@f0(theta = theta.temp[k,],
                                            args = c(f0.args$common.args,
                                            f0.args$pattern.args[[i]][[j]]))
                      }
                      
                      if (patternObj@ordered & length(hypotheses[[i]])>1)
                      {
                        zhat[,k,i] <- foo+family@deflate(theta = theta.temp[k,],f0.args$pattern.args[[i]])
                      } else {
                        zhat[,k,i] <- foo
                      }
                    }
                  }

                  #update posterior prob
                  zhat <- .Call("makeProbVectArr", zhat, p.temp, PACKAGE = "EBarrays")

                  #complete.lik<-zhat #value unimportant, just need the same dimension
                  #for (i in 1:n.gene)
                  #{
                  #  complete.lik[i,,]<-exp(zhat[i,,])*p.temp
                  #}
                  #
                  #zhat<-complete.lik/apply(complete.lik,1,sum)
                  
                  p.temp<-apply(zhat,c(2,3),mean)

                  if (independent)
                  {
                      p.pat<-apply(p.temp,1,sum)
                      p.clus<-apply(p.temp,2,sum)
                      p.temp<-outer(p.pat,p.clus,"*")
                  }
                          
                  ## M step
                  bic<-0
                  for (i in 1:n.cluster) {
                      zhatlogp <- sum(zhat[,i,] %*% log(p.temp[i,]))
                      max.comp.lik <-
                          optim(par = theta.temp[i,],
                                fn = complete.loglik,
                                method = c("L-BFGS-B"),
                                lower = family@lower.bound,
                                upper = family@upper.bound,
                                control = optim.control,
                                zhat = zhat[,i,],
                                zhatlogp = zhatlogp,
                                hypotheses = hypotheses,
                                f0 = family@f0,
                                f0.args = f0.args)
                    
                      theta.temp[i,] <- max.comp.lik$par
                      bic<-bic+max.comp.lik$value
                  }

                  if (verbose > 1) {
                      cat("\n\n  theta.temp: ")
                      print(theta.temp)
                  }
              }

              ## Finish up
              if (verbose > 1) {
                  cat("\n\n")
                  print("model parameter estimates:")
                  print(theta.temp)
                  print("Estimates of mixing proportions:")
                  print(t(p.temp))
              }

              if (verbose && exists("proc.time")) {
                  endTime <- proc.time()
                  cat(sprintf("\n\n Fit used %5.2f seconds user time\n",
                              endTime[1] - startTime[1]))
              }

              ic<- bic+(3*n.cluster+length(hypotheses)*n.cluster-1)*c(0.5*log(sum(all.pos)),1,log(log(sum(all.pos))))
              names(ic)<-c("BIC","AIC","HQ")
              
              new("ebarraysEMfit", 
                   family = family,
                   hypotheses = patternObj,
                   thetaEst = theta.temp,
                   probEst = p.temp,
                   ic = ic)
}

emfit2<-function(data,
                 family,
                 hypotheses,
                 cluster,
                 num.iter = 20,
                 verbose = getOption("verbose"),
                 optim.control = list(),
                 ...) {

              if (is(data, "ExpressionSet")) data<-exprs(data)
              
              family <- as(family, "ebarraysFamily")
              
              if (length(cluster)!=dim(data)[1]) {
                stop("Number of genes appeared in cluster does not match the data!")
              }

              thetaEst <- NULL
              probEst <- NULL
              
              unique.cluster<-unique(cluster)

              ic<-rep(0,3)
              for (i in unique.cluster) {
                fit <- emfit0(data,
                              family,
                              hypotheses,
                              num.iter,
                              verbose,
                              optim.control,
                              ...)

                ic<-ic+fit@ic
                
                thetaEst<-rbind(thetaEst,fit@thetaEst[1,])
                probEst<-rbind(probEst,fit@probEst[1,])
              }
              
              rownames(thetaEst)<-unique.cluster
              rownames(probEst)<-unique.cluster
              
              new("ebarraysEMfit", 
                   family = family,
                   hypotheses = hypotheses,
                   thetaEst = thetaEst,
                   probEst = probEst,
                   ic = ic)
}

emfit3<-function(data,
                 family,
                 hypotheses,
                 n.cluster,
                 criterion="BIC",
                 cluster.init = NULL,
                 independent = FALSE,
                 num.iter = 20,
                 verbose = getOption("verbose"),
                 optim.control = list(),
                 ...) {

             crit<-switch(EXPR=criterion,BIC=1,AIC=2,HQ=3)
             n.cluster<-sort(n.cluster)
             fit<-emfit1(data,family,hypotheses,n.cluster[1],cluster.init,independent, num.iter,verbose,optim.control)

             for (i in 2:length(n.cluster)) {
               cur.fit<-emfit1(data=data,family=family,hypotheses=hypotheses,n.cluster=n.cluster[i],cluster.init=cluster.init,independent=independent, num.iter=num.iter,verbose=verbose,optim.control=optim.control)

               if (cur.fit@ic[crit]<fit@ic[crit]) {
                 fit<-cur.fit
               }
             }

             return(fit)
           }

emfitmv<-          function(data,
                   family,
                   hypotheses,
                   groupid = NULL,
                   theta.init = NULL, p.init = NULL,
                   num.iter = 20,
                   verbose = getOption("verbose"),
                   trace = TRUE,
                   optim.control = list(),
                   ...) {

              ## Store current process time
              
              if (hypotheses@ordered)
                 stop("Family LNNMV does not support ordered hypotheses.")
                 
              if (exists("proc.time")) startTime <- proc.time()

              ## Define function to be optimized by optim:

              complete.loglik <-
                  function(theta,
                           zhat,
                           zhatlogp,
                           hypotheses,
                           f0, f0.args)
                  {
                      pred.dens <- zhat # values unimportant, just need the same dimensions

                      ## NOTE: f0 returns the log value of the predictive probability
                      ## density function.
                      for (i in seq(along = hypotheses))
                      {
                          foo <- 0 #numeric(nrow(zhat))
                          for (j in seq(along = hypotheses[[i]]))
                          {
                              foo <- 
                                  foo +
                                      f0(theta = theta,
                                         args = c(f0.args$common.args,
                                         f0.args$pattern.args[[i]][[j]]))
                          }
                          pred.dens[, i] <- foo
                      }
                      -(sum(zhat * pred.dens) + zhatlogp)
                  }


              ## Some preprocessing:
              verbose <- as.numeric(verbose)
              patternObj <- hypotheses
              hypotheses <- as(hypotheses@patterns, "list")

              if (verbose > 0) cat(paste("\n Checking for negative entries..."))
              all.pos <- apply(data, 1, function(x) !any(x <= 0))
              if (verbose > 0 && sum(all.pos) != length(all.pos)) {
                  cat(paste("\n\t", length(all.pos) - sum(all.pos),
                            "rows out of", length(all.pos),
                            "had at least one negative entry"))
                  cat("\n\t These rows will not be used in the EM fit\n")
              }


              if (verbose > 0)
                  cat(paste("\n Generating summary statistics for patterns.",
                            "\n This may take a few seconds...\n"))

              if (family=="LNNMV"){
                      f0.args <- 
                          family@f0.arglist(data[all.pos, ],
                                            hypotheses,groupid)
              }
              else
                f0.args <- 
                          family@f0.arglist(data[all.pos, ],
                                            hypotheses)
              
              

              theta.temp <-

                  if (is.null(theta.init))

                      ## note that this is expected in terms of the
                      ## parametrization used in actual optimization,
                      ## not necessarily the one expected by the user.

                      family@thetaInit(data[all.pos, ])

                  else

                      ## explicitly supplied initial estimates are in
                      ## the 'user' parametrization, which are
                      ## converted to the required parametrization by
                      ## the family@link function

                      family@link(theta.init)

              ## should we rm(data) now ? Is data copied in memory ?

              theta.mat <- matrix(NA, num.iter + 1, length(theta.temp))
              theta.mat[1,] <- theta.temp

              ## If missing, set initial probabilities as all equal
              if (is.null(p.init)) p.init <- rep(1, length(hypotheses))
              p.init <- p.init / sum(p.init)
              p.temp <- p.init
              p.mat <- matrix(NA, num.iter + 1, length(p.init))
              p.mat[1,] <- p.init

              ## Start Iterations
              if (verbose > 0)
                  cat(paste("\n Starting EM iterations (total", num.iter,
                            ").\n This may take a while\n\n"))

              notdone <- TRUE ## not used now, perhaps some convergence criteria
              iter <- 0
              zhat <- matrix(0, sum(all.pos), length(hypotheses))

              while ({iter <- iter + 1;  (iter <= num.iter && notdone)})
              {

                  if (verbose > 0)
                      cat(paste("\t Starting iteration", iter, "...\n"))
                  
                  ## E step
                  zhat[,] <- 0

                  ## NOTE: f0 returns the log value of the predictive
                  ## probability density function. Thus to get the log of the
                  ## product of the predictive probability densities, we add
                  ## (log(PA=1 * PA=2) = log(PA=1) + log(PA=2)).

                  for (i in seq(along = hypotheses))
                  {
                      foo <- 0
                      for (j in seq(along = hypotheses[[i]]))
                      {
                          foo <- 
                              foo + 
                                  family@f0(theta = theta.temp,
                                            args = c(f0.args$common.args,
                                            f0.args$pattern.args[[i]][[j]]))
                      }
                      zhat[, i] <- foo
                  }

                  zhat <- .Call("makeProbVect", zhat, p.temp, PACKAGE = "EBarrays")

                  p.temp <- colSums(zhat)
                  p.temp <- p.temp / sum(p.temp)
                  if (verbose > 1) {
                      cat("\n\n      p.temp: ")
                  }


                  ## M step
                  zhatlogp <- sum(zhat %*% log(p.temp))
                  theta.temp <-
                      optim(par = theta.temp,
                            fn = complete.loglik,
                            method = c("L-BFGS-B"),
                            lower = family@lower.bound,
                            upper = family@upper.bound,
                            control = optim.control,
                            zhat = zhat,
                            zhatlogp = sum(zhat %*% log(p.temp)),
                            hypotheses = hypotheses,
                            f0 = family@f0,
                            f0.args = f0.args)$par

                  theta.mat[iter + 1,] <- theta.temp
                  p.mat[iter + 1,] <- p.temp

                  if (verbose > 1) {
                      cat("\n\n  theta.temp: ")
                      print(theta.temp)
                  }
              }

              ################################
              rm(zhat)
              ################################

              ## Finish up

              colnames(theta.mat) <- names(theta.temp) <-
                  c(paste("theta", 1:(length(theta.temp)), sep=""))
              colnames(p.mat) <- names(p.temp) <-
                  c(paste("p", 1:(length(p.temp)), sep=""))

              if (verbose > 1) {
                  cat("\n\n")
                  print("model parameter estimates:")
                  print(theta.temp)
                  print("Estimates of mixing proportions:")
                  print(p.temp)
              }

              if (verbose && exists("proc.time")) {
                  endTime <- proc.time()
                  cat(sprintf("\n\n Fit used %5.2f seconds user time\n",
                              endTime[1] - startTime[1]))
              }

              new("ebarraysEMfit", 
                   family = family,
                   hypotheses = patternObj,
                   thetaEst = rbind(theta.temp,NULL),
                   probEst = rbind(p.temp,NULL),
                   ic = 0)
}


#posterior probability calculation
setGeneric("postprob",
           function(fit,data, ...) standardGeneric("postprob"))

setMethod("postprob",
          signature(fit = "ebarraysEMfit",
                    data = "matrix"),

          function(fit,
                   data,...) {

    data[data <= 0] <-min(data[data > 0])
    patternObj<-fit@hypotheses
    hypotheses <- patternObj@patterns
    family<-fit@family

    n.cluster<-dim(fit@thetaEst)[1]
    
    f0.args <-family@f0.arglist(data,hypotheses,...)
    zhat<-array(0,dim=c(dim(data)[1], n.cluster,length(hypotheses)))
                  
    ## NOTE: f0 returns the log value of the predictive
    ## probability density function. Thus to get the log of the
    ## product of the predictive probability densities, we add
    ## (log(PA=1 * PA=2) = log(PA=1) + log(PA=2)).
    for (k in 1:n.cluster)
    {
        for (i in seq(along = hypotheses))
        {
            foo <- 0
            for (j in seq(along = hypotheses[[i]]))
            {
                foo <- 
                    foo + 
                        family@f0(theta = fit@thetaEst[k,],
                                  args = c(f0.args$common.args,
                                  f0.args$pattern.args[[i]][[j]]))
            }
            
            if (patternObj@ordered & length(hypotheses[[i]])>1)
            {
                zhat[,k,i] <- foo+family@deflate(theta = fit@thetaEst[k,],f0.args$pattern.args[[i]])
            } else {
                zhat[,k,i] <- foo
            }
        }
    }

    #update posterior prob
    zhat <- .Call("makeProbVectArr", zhat, fit@probEst, PACKAGE = "EBarrays")

    #complete.lik<-zhat #value unimportant, just need the same dimension
    #for (i in 1:dim(data)[1])
    #{
    #    complete.lik[i,,]<-exp(zhat[i,,])*fit@probEst
    #}    
    #post.prob<-complete.lik/apply(complete.lik,1,sum)

    rownames(zhat)<-rownames(data)

    return(list(joint=zhat,cluster=apply(zhat,c(1,2),sum),pattern=apply(zhat,c(1,3),sum)))
})
