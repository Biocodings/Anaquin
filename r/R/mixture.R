#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

loadFusMix <- function(data, file='A')
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

loadVarMix <- function(data, file='A')
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

loadTransMix <- function(file='A')
{
    if (is.null(file))
    {
        return(NULL)
    }

    stopifnot(identical(file, 'A') || identical(file, 'B') || identical(file, 'F') || identical(file, 'AB'))

    if (identical(file, 'A'))  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR002.v013.csv') }
    if (identical(file, 'B'))  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR003.v013.csv') }
    if (identical(file, 'F'))  { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR005.v013.csv') }
    if (identical(file, 'AB')) { file <- url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv') }
    
    mix <- read.csv(file, row.names=1, sep='\t')
    
    stopifnot(ncol(mix) == 2 || ncol(mix) == 3)
    isAB <- ncol(mix) == 3

    if (isAB) { colnames(mix) <- c('length', 'A', 'B') } else { colnames(mix) <- c('length', 'A') }

    #
    #           length     A           B
    # R1_11_1    703  161.13281    5.035400
    # R1_11_2    785   80.56641   10.070801
    #
    
    if (isAB)
    {
        mix$fold    <- mix[[3]] / mix[[2]]
        mix$logFold <- log(mix$fold)
    }
    
    # Eg: R1_1 for R1_1_1 and R1_1_2    
    mix$GeneID <- isoformsToGenes(row.names(mix))
    
    # Genes that are defined in the mixture
    gIDs <- unique(mix$GeneID)

    if (isAB)
    {
        gd <- data.frame(row.names=gIDs,
                         A=rep(0, length(gIDs)),
                         B=rep(0, length(gIDs)),
                         fold=rep(0, length(gIDs)),
                         logFold=rep(0, length(gIDs)))
    }
    else
    {
        gd <- data.frame(row.names=gIDs, A=rep(0, length(gIDs)))
    }
        
    for (id in gIDs)
    {
        seqs <- mix[mix$GeneID == id,]

        #
        #      ID     Length     Mix A      Mix B
        #    -------------------------------------
        #     R1_11    703     161.13281    5.0354
        #
        
        gd[id,]$A <- sum(seqs[,2])

        if (isAB)
        {
            gd[id,]$B <- sum(seqs[,3])
            gd[id,]$fold <- gd[id,]$B / gd[id,]$A
            gd[id,]$logFold <- log2(gd[id,]$fold)
        }
    }
    
    i <- data.frame(mix)
    r <- list('isoforms'=i[with(i, order(row.names(i))),], 'genes'=gd[with(gd, order(row.names(gd))),])

    class(r) <- c("TransMixture")
    return (r)
}