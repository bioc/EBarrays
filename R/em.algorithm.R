
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


## Class for result of emfit : ebarraysEmFit

setClass("ebarraysEmFit",
         representation(family = "ebarraysFamily",
                        hypotheses = "ebarraysPatterns",
                        thetaEst = "numeric",
                        probEst = "numeric",
                        thetaTrace = "matrix",
                        probTrace = "matrix"))


setMethod("show", "ebarraysEmFit",

          function(object) {

              cat(paste("\n EB model fit",
                        "\n\t Family:", object@family,
                        "(", object@family@description, ")"))
              cat("\n\n Model parameter estimates:\n\n")
              print(object@family@invlink(object@thetaEst))
              cat("\n Estimated mixing proportions:\n\n")
              print(object@probEst)
              cat("\n\n Additional slots: @hypotheses, @thetaTrace, @probTrace\n\n")

              invisible(object)
          })


setGeneric("emfit",
           function(data, family, hypotheses, ...) standardGeneric("emfit"))


setMethod("emfit",
          signature(data = "exprSet",
                    family = "character",
                    hypotheses = "ebarraysPatterns"),

          function(data,
                   family,
                   hypotheses, ...) {

              data <- exprs(data)
              family <- as(family, "ebarraysFamily")
              callGeneric(data = data,
                          family = family,
                          hypotheses = hypotheses, ...)

          })



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
          signature(data = "exprSet",
                    family = "ebarraysFamily",
                    hypotheses = "ebarraysPatterns"),

          function(data,
                   family,
                   hypotheses, ...) {

              data <- exprs(data)
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
                   theta.init = NULL, p.init = NULL,
                   num.iter = 20,
                   verbose = getOption("verbose"),
                   trace = TRUE,
                   optim.control = list(),
                   ...) {

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
              hypotheses <- as(hypotheses, "list")

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

              if (trace)
                  new("ebarraysEmFit", 
                      family = family,
                      hypotheses = patternObj,
                      thetaEst = theta.temp,
                      probEst = p.temp,
                      thetaTrace = theta.mat,
                      probTrace = p.mat)
              else
                  new("ebarraysEmFit", 
                      family = family,
                      hypotheses = patternObj,
                      thetaEst = theta.temp,
                      probEst = p.temp)
          })


