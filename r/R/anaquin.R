#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.createData <- function(x, keys)
{
    # Do we have sequin names?
    if (is.null(x$names))
    {
        return (NULL)
    }
    
    data <- data.frame(row.names=x$names)
    
    for (key in keys)
    {
        if (!is.null(x[[key]]))
        {
            data[[key]] <- x[[key]]
        }
    }

    return (data)    
}

.createMixture <- function(x)
{
    if (!is.null(x))
    {
        #
        # There're several possibilities:
        #
        #     1. local file (eg: /home/tedwon/data/mixture.csv)
        #     2. URL (eg: https://s3.amazonaws.com/anaquin/mixtures/MTR002.v013.csv)
        #     3. Data frame (or matrix)
        #
    
        if (is(x, 'url'))        { return (read.csv(x)) }
        if (is(x, 'data.frame')) { return (x) }
        if (is(x, 'matrix'))     { return (data.frame(x)) }
        
        if (identical(x, 'A'))  { return('A')  }
        if (identical(x, 'B'))  { return('B')  }
        if (identical(x, 'F'))  { return('F')  }
        if (identical(x, 'AB')) { return('AB') }
        
        if (is(x, 'character'))  { return (read.csv(x)) }
    }
    
    return (NULL)
}

CreateDataForAnaquin <- function(...)
{
    x <- list(...)

    # Metrics Anaquin supports
    keys <- c('sd', 'pval', 'qval', 'ratio', 'mean', 'expected', 'input', 'measured', 'sensitivity', 'label', 'score')

    r <- list('seqs'=.createData(x, keys))
    
    if (is.null(r$seqs))
    {
        stop('Please specify sequin names, for example, CreateDataForAnaquin(names=...)')
    }
    
    class(r) <- 'Anaquin'
    return (r)
}

########################################################
#                                                      #
#                  RnaQuin functions                   #
#                                                      #
########################################################

#
# Eg: R1_1_1 to R1_1. Works for TransQuin sequins.
#
isoformsToGenes <- function(ids)
{
    ids <- strsplit(ids, '_')
    
    f <- function(x)
    {
        paste(x[1:length(x)-1], collapse='_')
    }
    
    unlist(lapply(ids, f))
}

#
# This function takes a TransQuin data set and a list of gene IDs. It gives the expected expression for
# each of the gene requested.
#

expectForGenes <- function(data, gIDs)
{
    stopifnot(class(data) == 'TransQuin')
    data <- data$seqs
        
    r <- data.frame(expected=rep(NA, length(gIDs)), efrac=rep(NA, length(gIDs)), row.names=gIDs)
    
    for (gID in row.names(r))
    {
        x <- data[grep(gID, row.names(data)),]
        
        stopifnot(nrow(x) >= 1)
        
        minor <- x[which.min(x$expected),]$expected
        major <- x[which.max(x$expected),]$expected
        
        r[row.names(r) == gID,]$efrac <- (minor / major)
        r[row.names(r) == gID,]$expected <- sum(x$expected)
    }
    
    return (r)
}

########################################################
#                                                      #
#                VarQuin functions                     #
#                                                      #
########################################################

#
# Collapse VarQuin sequins to sequin genes. For example, 'D_2_8_R' and 'D_2_8_V' to just 'D_2_8'.
#

VarQuin.genes <- function(seqs)
{
    seqs <- gsub('_R', '', seqs)
    seqs <- gsub('_V', '', seqs)
    
    return (unique(seqs))
}

alleleFreq <- function(data, ids)
{
    stopifnot(class(data) == 'VarQuin' | class(data) == 'VarMixture')
    
    if (class(data) == 'VarQuin')
    {
        data <- data$mix
    }
    
    freq <- data$allFreq
    freq <- freq[row.names(freq) %in% ids,]
    
    return (freq)
}

########################################################
#                                                      #
#                FusQuin functions                     #
#                                                      #
########################################################

normalFusion <- function(data)
{
    stopifnot(class(data) == 'FusQuin' | class(data) == 'FusMixture')
    
    if (class(data) == 'FusQuin')
    {
        data <- data$mix
    }

    require(dplyr)

    x <- mix$seqs[mix$seqs$label!='T',]
    x$names <- row.names(x)
    
    setDT(x)[, Grp:= .GRP, label]
    x <- dcast(x[, N:= 1:.N, label], N~Grp, value.var=c('names', 'abund'), sep='')[,N:= NULL][]
    x$abund <- x$abund1 / x$abund2
        
    return (data.frame(x))
}

normalToFusion <- function(names)
{
    return (cat('F', substr(names, 2, nchar(names)), sep=''))
}