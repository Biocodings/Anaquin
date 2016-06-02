#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.createData <- function(x, keys)
{
    if (is.null(x$seqs))
    {
        return (NULL)
    }
    
    data <- data.frame(row.names=x$seqs)
    
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

TransQuin <- function(...)
{
    x <- list(...)
    
    #
    # seqs      Unique sequin names
    # input:    Expected abundance (attomol/ul etc) 
    # measured: Measured abundance
    # mean:     Replicate means
    # pval:     Probability under the null hypothesis
    # qval:     Adjusted probability under the null hypothesis
    #
    keys <- c('pval', 'qval', 'mean', 'expected', 'input', 'measured')

    r <- list('seqs'=.createData(x, keys))
    
    if (is.null(r$seqs))
    {
        r$mix <- loadTransMix(.createMixture(x$mix))
        
        if (!is.null(r$mix))
        {
            # This is differential AB...
            if (ncol(r$mix$isoforms) == 6)
            {
                r$seqs <- data.frame(row.names=row.names(r$mix$isoforms), expected=r$mix$isoforms$fold)
            }
        }
        else
        {
            stop('No sequin found. Please check and try again.') 
        }
    }

    stopifnot(!is.null(r$seqs))
    
    class(r) <- 'TransQuin'
    return (r)
}

MetaQuin <- function(...)
{
    x <- list(...)
    
    keys <- c('mix', 'seqs', 'input', 'expected', 'measured')
    data <- .createData(x, keys)

    r <- list('seqs'=data, mix=.createMixture(x$mix))
    
    if (!is.null(x[['bedgr']])) { r$bedgr <- x$bedgr }
    if (!is.null(x[['annot']])) { r$annot <- x$annot }
    
    stopifnot(!is.null(r$seqs))
    
    class(r) <- 'MetaQuin'
    return (r)
}

VarQuin <- function(...)
{
    x <- list(...)
    
    #
    # seqs     Unique sequin names
    # label    Sequin classification
    # rRead    Number of reads supporting the reference allele
    # vRead    Number of reads supporting the alternative allele
    # type     Variant type
    # ratio    Input fold-change
    # input    Input concentration (eg: allele frequency)
    # measured Measured expression (eg: allele frequency)
    # bedgr    Bedgraph file name
    # annot    Annotation file name
    #
    
    keys <- c('label', 'pval', 'rRead', 'vRead', 'type', 'ratio', 'input', 'measured', 'bedgr', 'annot')
    data <- .createData(x, keys)
    
    r <- list('seqs'=data, mix=.createMixture(x$mix))
    
    if (!is.null(x[['bedgr']])) { r$bedgr <- x$bedgr }
    if (!is.null(x[['annot']])) { r$annot <- x$annot }

    # TODO: This should be activated, but won't work for VarCoverage_density.R (we should read the sequins from the file directly)
    #stopifnot(!is.null(r$seqs))
    
    class(r) <- 'VarQuin'
    return (r)
}

FusQuin <- function(...)
{
    x <- list(...)
    
    #
    # seqs      Unique sequin names
    # input:    expected abundance (attomol/ul etc) 
    # measured: measured abundance
    # label:    sequin classification
    #

    keys <- c('input', 'measured', 'label')
    r <- list('seqs'=.createData(x, keys))
    
    if (is.null(r$seqs))
    {
        r$mix <- loadFusMix(.createMixture(x$mix))
        
        if (is.null(r$mix))
        {
            stop('No mixture found. Please check and try again.') 
        }
    }

    stopifnot(!is.null(r$seqs))
    
    class(r) <- 'FusQuin'
    return (r)
}

StructQuin <- function(...)
{
    x <- list(...)
    
    keys <- c()
    data <- .createData(x, keys)
    
    r <- list('seqs'=data, mix=.createMixture(x$mix))
    
    stopifnot(!is.null(r$seqs))
    
    class(r) <- 'StructQuin'
    return (r)    
}

LadQuin <- function(...)
{
    x <- list(...)
    
    #
    # seqs     Unique sequin names
    # input    Input abundance (eg: attomol/ul)
    # measured Measured abundance
    # mix      Mixture details
    #
    
    keys <- c('label', 'input', 'measured')
    data <- .createData(x, keys)
    
    r <- list('seqs'=data, mix=.createMixture(x$mix))

    stopifnot(!is.null(r$seqs))

    class(r) <- 'LadQuin'
    return (r)
}

########################################################
#                                                      #
#                TransQuin functions                   #
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