#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

.isoformsToGenes <- function(trans)
{
    trans <- as.character(trans)
    genes <- substr(as.character(trans), 1, nchar(trans)-2)
    genes
}

loadGene <- function(id, mix=loadMixture())
{
    r <- mix$genes[row.names(mix$genes)==id,]
    r
}

#
# Returns a list of exon bins that can be used as negative control for normalization (ie: FC==1)
#

negativeExonBins <- function(m = loadMixture())
{
    exons <- m$exons[m$exons$fold==1,]

    # Make sure the corresponding logFolds are correct
    stopifnot(nrow(exons[exons$logFold != 0,]) == 0)

    #
    # For example, the names of bins can be accessed by row.names(exons)
    #
    
    exons
}

#
# Returns a list of isoforms that can be used as negative control for normalization (ie: FC==1)
#

negativeIsoforms <- function(m = loadMixture())
{
    i <- m$isoforms[m$isoforms$fold==1,]
    
    # Make sure the corresponding logFolds are correct
    stopifnot(nrow(i[i$logFold != 0,]) == 0)
    
    #
    # For example, the names of bins can be accessed by row.names(exons)
    #
    
    i
}

expExonBins <- function(m = loadMixture())
{
    r <- negativeExonBins()
    r <- r[c('R2_71:E001', 'R2_71:E002', 'R2_71:E003', 'R2_71:E004', 'R2_71:E005', 'R2_71:E006', 'R2_71:E007', 'R2_71:E008'),]
    r
}

#
# Load the mixture into an R object that can be used in other Anaquin functions
#

loadMixture <- function(mix=NULL, exons=NULL)
{
    if (is.null(mix))
    {
        mix <- read.csv(url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv'), row.names=1, sep='\t')
        colnames(mix) <- c('Length', 'Mix.A', 'Mix.B')
    }

    if (is.null(exons))
    {
        #exons <- read.csv(url('https://s3.amazonaws.com/anaquin/mixtures/MTR006.v013.csv'), row.names=1)
        #exons$logFold <- as.numeric(as.character(exons$logFold))
    }
    
    mix$fold    <- mix$Mix.B / mix$Mix.A
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

    r <- list('isoforms'=i, 'genes'=g[with(g, order(row.names(g))),], 'exons'=exons)
    class(r) <- c("Mixture")
    
    r$genes$logFold <- round(r$genes$logFold)
    r
}
