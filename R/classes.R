

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


## Define miscellaneous classes 


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
                        upper.bound = "numeric"))

setAs("character", "ebarraysFamily",

      function(from) {

          if (from == "GG") 
              eb.createFamilyGG()
          else if (from == "LNN")
              eb.createFamilyLNN()
          else stop(paste("\n The only families recognized by name are GG (Gamma-Gamma)",
                          "\n and LNN (Lognormal-Normal). Other families need to be",
                          "\n specified as objects of class 'ebarraysFamily'"))
          
      })


setClass("ebarraysPatterns",
         representation("list"))


setMethod("show", "ebarraysFamily",

          function(object) {
              cat(paste("\n Family:", object, "(", object@description, ")",
                        "\n\nAdditional slots:\n"))
              show(getSlots(class(object)))
              invisible(object)
          })


setMethod("show", "ebarraysPatterns",

          function(object) {
              cat(paste("\n Collection of", length(object), "patterns"))
              for (i in seq(along = object))
              {
                  pat = object[[i]]
                  cat("\n\n    Pattern", i, "has", length(pat),
                      if (length(pat) == 1) "group" else "groups")
                  for (j in seq(along = pat))
                  {
                      cat(sprintf("\n       Group %2d: ", j))
                      cat(paste(pat[[j]], collapse = " "))
                  }
              }
              cat("\n")
              invisible(object)
          })


