#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

#library("Rsamtools")
#library("GenomicFeatures")
#library("GenomicAlignments")

.isoformsToGenes <- function(trans)
{
    trans <- as.character(trans)
    genes <- substr(as.character(trans), 1, nchar(trans)-2)
    genes
}

sequin <- function(id, mix=loadMixture())
{
    r <- mix$genes[mix$genes==id,]
    r
}

fold <- function(d)
{
    r <- d$Fold
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
        mix <- read.csv(url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv'), sep='\t')
        colnames(mix) <- c('ID', 'Length', 'Mix.A', 'Mix.B')
    }

    if (is.null(exons))
    {
        exons <- read.csv('/Users/tedwong/Desktop/exons.csv', row.names=1)
        exons$logFold <- as.numeric(as.character(exons$logFold))
    }

    # Eg: R1_1 for R1_1_1 and R1_1_2    
    mix$GeneID <- .isoformsToGenes(mix$ID)
    
    # Genes that are defined in the mixture
    geneIDs <- unique(mix$GeneID)
    
    g <- data.frame(ID=geneIDs,
                    A=rep(0, length(geneIDs)),
                    B=rep(0, length(geneIDs)),
                    Fold=rep(0, length(geneIDs)),
                    LogFold=rep(0, length(geneIDs)),
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
        
        g[g$ID == id,]$A <- sum(seqs[,3])
        g[g$ID == id,]$B <- sum(seqs[,4])
        
        #
        # Calculate the expected fold change (B/A)
        #
        
        g[g$ID == id,]$Fold    <- g[g$ID == id,]$B / g[g$ID == id,]$A
        g[g$ID == id,]$LogFold <- log2(g[g$ID == id,]$Fold)
        
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
        
        g[g$ID == id,]$type <- type
    }
    
    # Prefer not to have it as factor variable
    g$ID <- as.character(g$ID)
    
    # Sort by ID so that it'll more easier interpreted
    g <- g[with(g, order(ID)),]
    
    i <- data.frame(mix)
    i <- i[with(i, order(ID)),]

    r <- list('isoforms'=i, 'genes'=g, 'exons'=exons)
    class(r) <- c("Mixture")
    r
}
