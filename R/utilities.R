
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


createExprSet <-
    function(datafile, repfile)
{
    ## Intensities
    data.temp <-
        read.table(datafile, header = TRUE,
                   row.names = NULL, comment.char = "")
    int.data <- data.temp[,-1]
    dimnames(int.data)[[1]] <- data.temp[,1]
    dimnames(int.data)[[2]] <- dimnames(data.temp)[[2]][-1]

    ## Phenodata
    replicates <- data.frame(replicates = factor(scan(repfile, quiet = TRUE)))

    ## initially planned to create exprSet, but let's skip that till
    ## we decide whether we want to depend on bioconductor. Let this
    ## be just a matrix for now, where the repfile is essentially
    ## ignored. Not really used anywhere anyway, so not a big deal.

    ## Create and return object of class exprSet 
    new("exprSet", exprs = as.matrix(int.data),
        phenoData = new("phenoData", pData = replicates,
        varLabels = list(replicates = "Replicates")))

    ## back to exprSet
    ##as.matrix(int.data)
}




ebPatterns <-
    function(x)
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
    new("ebarraysPatterns", patterns)
}


