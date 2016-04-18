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
        if (is(x, 'character'))  { return (read.csv(x)) }
        if (is(x, 'data.frame')) { return (x) }
        if (is(x, 'matrix'))     { return (data.frame(x)) }
    }
    
    return (NULL)
}

TransQuin <- function(...)
{
    x <- list(...)
    
    keys <- c('label', 'pval', 'qval', 'mean', 'se', 'ratio', 'expected', 'measured')
    data <- .createData(x, keys)

    r <- list('seqs'=data, mix=.createMixture(x$mix))
    class(r) <- 'TransQuin'
    
    return (r)
}

VarQuin <- function(...)
{
    x <- list(...)
    
    keys <- c('label', 'pval', 'rRead', 'vRead', 'type', 'ratio', 'expected', 'measured', 'bedgr', 'annot')
    data <- .createData(x, keys)
    
    r <- list('seqs'=data, mix=.createMixture(x$mix))
    class(r) <- 'VarQuin'
    
    if (!is.null(x[['bedgr']])) { r$bedgr <- x$bedgr }
    if (!is.null(x[['annot']])) { r$annot <- x$annot }

    return (r)
}

FusQuin <- function(...)
{
    x <- list(...)
    
    keys <- c('label', 'pos1', 'pos2')
    data <- .createData(x, keys)

    r <- list('seqs'=data, mix=.createMixture(x$mix))
    class(r) <- 'FusQuin'
    
    return (r)
}

LadQuin <- function(...)
{
    x <- list(...)
    
    keys <- c('label', 'elfc', 'lfc', 'pval', 'abund', 'type', 'expected', 'measured', 'aligned')
    data <- .createData(x, keys)
    
    r <- list('seqs'=data, mix=.createMixture(x$mix))
    class(r) <- 'LadQuin'
    
    return (r)
}

########################################################
#                                                      #
#                 Online Mixtures                      #
#                                                      #
########################################################

#
# Returns a URL of mixture for TransQuin. This saves the trouble of downloading it explicitly.
#

transMixURL <- function()
{
    
}

########################################################
#                                                      #
#                Accessor functions                    #
#                                                      #
########################################################

#
# Filtering functionality The following are supported:
#
#   'seqs'
#   'sequins'
#   'endos'
#
#filter <- function(data, name)
#{
#    stopifnot(class(data) == 'TransQuin' |
#                  class(data) == 'VarQuin'   |
#                  class(data) == 'MetaQuin')
#    
#    stopifnot(name == 'seqs'    |
#                  name == 'sequins' |
#                  name == 'endos')
#    
#    if (name == 'seqs' | name == 'sequins')
#    {
#        if (!is.null(data$seqs$elf))
#        {
#            return (data$seqs[!is.na(data$seqs$elf),])
#        }
#        else
#       {
#            return (data$seqs) # TODO: Fix this
#        }
#    }
#    else if (name == 'endos')
#    {
#        return (data$seqs[is.na(data$seqs$elf),])        
#    }
#}

#
# Normalize the counts and calculate the base mean. Base mean is the average of the normalized count values, taken
# over all samples.
#

#baseMean <- function(data)
#{
#    stopifnot(class(data) == 'TransQuin')
#    
#    if (is.null(data$seqs$baseMean))
#    {
#        # TODO: Implement me
#    }
#    
#    return (data$seqs$baseMean)
#}

seqs <- function(data)
{
    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin'  |
              class(data) == 'LadQuin')

    return (row.names(data$seqs))
}

.isoformsToGenes <- function(trans)
{
    trans <- as.character(trans)
    genes <- substr(as.character(trans), 1, nchar(trans)-2)
    genes
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

########################################################
#                                                      #
#                 Mixture Loading                      #
#                                                      #
########################################################

loadMixture.FusQuin <- function(file='/Users/tedwong/Sources/QA/data/FusQuin/MFU007.v013.csv')
{
    data <- read.csv(file, row.names=1, sep='\t')
    colnames(data) <- c('length', 'abund')

    data$label <- 'T'
    data[grep("FG", rownames(data)),]$label <- 'F'
    data[grep("NG", rownames(data)),]$label <- 'N'

    r <- list(seqs=data)
    class(r) <- c("FusMixture")
    
    return (r)
}

loadMixture.VarQuin <- function(file='A')
{
    if (file == 'A')  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MVA011.v013.csv') }
    if (file == 'F')  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MVA012.v013.csv') }

    mix <- read.csv(file, row.names=1, sep='\t')
    colnames(mix) <- c('length', 'A')
    
    #          length    A
    # D_1_1_R   1107  5002.451
    # D_1_2_R   1122  6669.968
    # D_1_3_R    986  8003.922    
    #

    seqs <- row.names(mix)
    seqs <- gsub('_R', '', seqs)
    seqs <- gsub('_V', '', seqs)
    seqs <- unique(seqs)
    
    # Table of expected allele frequencies
    allFreq <- data.frame(row.names=seqs, freq=rep(NA,length(seqs)))

    for (seq in seqs)
    {
        ref <- paste(seq, '_R', sep='')
        var <- paste(seq, '_V', sep='')        

        ref <- mix[row.names(mix) == ref,][[2]]
        var <- mix[row.names(mix) == var,][[2]]        
        
        allFreq[row.names(allFreq) == seq,] <- var / (ref + var)
    }
    
    r <- list(seqs=mix, allFreq=allFreq)
    class(r) <- c("VarMixture")

    return (r)
}

loadMixture.TransQuin <- function(file='A')
{
    if (file == 'A')  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv') }
    if (file == 'B')  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv') }
    if (file == 'F')  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv') }
    if (file == 'AB') { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv') }
    
    mix <- read.csv(file, row.names=1, sep='\t')
    colnames(mix) <- c('length', 'A', 'B')
    
    #
    #           length     A           B
    # R1_11_1    703  161.13281    5.035400
    # R1_11_2    785   80.56641   10.070801
    #

    mix$fold    <- mix[[3]] / mix[[2]]
    mix$logFold <- log(mix$fold)
    
    # Eg: R1_1 for R1_1_1 and R1_1_2    
    mix$GeneID <- .isoformsToGenes(row.names(mix))
    
    # Genes that are defined in the mixture
    geneIDs <- unique(mix$GeneID)
    
    g <- data.frame(row.names=geneIDs,
                    A=rep(0, length(geneIDs)),
                    B=rep(0, length(geneIDs)),
                    fold=rep(0, length(geneIDs)),
                    logFold=rep(0, length(geneIDs)),
                    type=rep(0, length(geneIDs)))
    
    #
    # Calculate the expected log-fold between mixture A and B
    #
    
    for (id in geneIDs)
    {
        seqs <- mix[mix$GeneID == id,]
        
        #
        # Calculate the expected abundance. We assume the following format:
        #
        #      ID     Length     Mix A      Mix B
        #    -------------------------------------
        #     R1_11    703     161.13281    5.0354
        #
        # 
        # We shouldn't assume anything for the column names.
        #
        
        g[id,]$A <- sum(seqs[,2])
        g[id,]$B <- sum(seqs[,3])
        
        #
        # Calculate the expected fold change (B/A)
        #
        
        g[id,]$fold    <- g[id,]$B / g[id,]$A
        g[id,]$logFold <- log2(g[id,]$fold)
        
        #
        # For each gene, we will classify it as:
        #
        #    'E' - equal concentration across mixtures
        #    'A' - concentration in the first mixture is higher
        #    'B' - concentration in the second mixture is higher
        #
        
        # We'll need that because we might have variable number of isoforms
        type = ''
        
        for (i in seqs$ID)
        {
            seq <- seqs[seqs$ID==i,]
            
            if (seq$Mix.A > seq$Mix.B)      { type <- paste(type, 'A', sep='') }
            else if (seq$Mix.B > seq$Mix.A) { type <- paste(type, 'B', sep='') }
            else                            { type <- paste(type, 'E', sep='') }
        }
        
        g[id,]$type <- type
    }
    
    i <- data.frame(mix)
    i <- i[with(i, order(row.names(i))),]
    
    r <- list('isoforms'=i, 'genes'=g[with(g, order(row.names(g))),])
    class(r) <- c("TransMixture")

    r$genes$logFold <- round(r$genes$logFold)
    r
}