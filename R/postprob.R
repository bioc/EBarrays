
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



setGeneric("postprob",
           function(fit, data, ...) standardGeneric("postprob"))



setClass("ebarraysPostProb",
         representation("matrix",
                        hypotheses = "ANY"))

setMethod("show", "ebarraysPostProb",

          function(object) {

              cat("Hypotheses:\n")
              show(object@hypotheses)

              cat(sprintf("\n\n Posterior Probabilities of %d hypotheses for %d genes:\n (summarized)\n\n",
                          ncol(object), nrow(object)))
              print(summary(as(object, "matrix")))

              invisible(object)
          })


##setMethod("[",
##          signature(x = "ebarraysPostProb",
##                    i = "ANY", j = "ANY", drop = "ANY"),
##          function(x, i, j, drop = FALSE, ...) {
##              x@.Data <- x@.Data[i, j, drop = drop, ...]
##              x
##          })



## setMethod("postprob",
##           signature(fit = "ebarraysSpFit", data = "missing"),

##           function(fit, data) {

##               new("ebarraysPostProb",
##                   fit@post,
##                   hypotheses = fit@hypotheses)
##           })






setMethod("postprob",
          signature(fit = "ebarraysEmFit", data = "exprSet"),

          function(fit, data) {

              data <- exprs(data)
              callGeneric(fit = fit, data = data)

          })



## add another (simple) method for the semiparametric fit later

setMethod("postprob",
          signature(fit = "ebarraysEmFit", data = "matrix"),

          function(fit, data) {

              data[data <= 0] <-
                  min(data[data > 0])
              hypotheses <- fit@hypotheses
              f0.args <-
                  fit@family@f0.arglist(data,
                                        hypotheses)

              pred.dens <- matrix(0, nrow(data), length(hypotheses))
              postprob <- matrix(NA, nrow(data), length(hypotheses),
                                 dimnames = 
                                 list(rownames(data),
                                      paste("P", seq(along = hypotheses), sep="")))
              
              for (i in seq(along = hypotheses))
              {
                  for (j in seq(along = hypotheses[[i]]))
                  {
                      pred.dens[, i] <-
                          pred.dens[, i] +
                              fit@family@f0.pp(theta = fit@thetaEst,
                                               args =
                                               c(f0.args$common.args,
                                                 f0.args$pattern.args[[i]][[j]]))
                  }
              }


              ## Is it worth moving this to C ?
              for (i in seq(along = hypotheses))
              {
                  out1.denom <- colSums(fit@probEst * t(exp(pred.dens - pred.dens[, i])))
                  postprob[, i] <- fit@probEst[i] / out1.denom
              }

              new("ebarraysPostProb", postprob, hypotheses = hypotheses)
          })



# setMethod("plot",
#           signature(x = "ebarraysPostProb"),
#           function(x, data = NULL, cutoff = .9, which = TRUE) {

#               if (is(data, "exprSet")) {
#                   print("Will eventually display patterns")
#                   show(x@patterns)
                  
#                   cat(sprintf("\n\n Posterior Probabilities of %d hypotheses for %d genes:\n (summarized)\n\n",
#                               ncol(x), nrow(x)))
#                   print(summary(as(x, "matrix")))
#               }
#               else
#                   show(x)

#               invisible(x)
#           })

