### Date: 20 March 2006
### Title: Empirical Bayes Analysis of Gene Expression Data
### Author: Ming Yuan <myuan@isye.gatech.edu>
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.

## Define classes for hypotheses

#Patterns for a single time point
setClass("ebarraysPatterns",
         representation(patterns = "list",
                        ordered = "logical"))

ebPatterns <-
    function(x,ordered=FALSE)
    ## x can be a character vector (of length > 2), a connection or a
    ## file name
{
    if (!(is.character(x) && length(x) > 1))
        x <- readLines(x)
    patterns <- vector("list", length(x))
    len <- FALSE
    for(i in seq(along = patterns))
    {
        pat <- as.numeric(strsplit(x[i], "\\W+")[[1]])
        if (is.logical(len))
        {
            len <- length(pat)
            if (len == 0) stop("Pattern has length 0")
        }
        if (length(pat) != len || any(is.na(pat)))
        {
            print(pat)
            stop("Invalid pattern")
        }
        vals <- sort(unique(pat[pat > 0]))
        patterns[[i]] <- vector("list", length(vals))
        for (j in seq(along = vals))
        {
            patterns[[i]][[j]] <- (1:len)[pat == vals[j]]
        }
    }
    new("ebarraysPatterns", patterns=patterns, ordered=ordered)
}

setMethod("show", "ebarraysPatterns",
          function(object) {
              cat(paste("\n Collection of", length(object@patterns), "patterns"))
              n.level<-0
              for (i in seq(along = object@patterns))
              {
                  pat = object@patterns[[i]]
                  n.level = max(n.level,length(pat))
                  
                  cat("\n\n    Pattern", i, "has", length(pat),
                      if (length(pat) == 1) "group" else "groups")
                  for (j in seq(along = pat))
                  {
                      cat(sprintf("\n       Group %2d: ", j))
                      cat(paste(pat[[j]], collapse = " "))
                  }
              }
              cat("\n")

              if (object@ordered)
              {
                  cat(paste("\n\n Patterns are ordered: Group 1 < ... < Group", n.level))
                  cat("\n")
              }
              
              invisible(object)
          })
