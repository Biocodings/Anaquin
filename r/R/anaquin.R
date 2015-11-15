#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")

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
# Load the mixture into an R object that can be used in other Anaquin functions
#

loadMixture <- function(mix)
{
    if (!hasArg(mix))
    {
        mix <- read.csv(url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv'), sep='\t')
    }

    # Eg: R1_1 for R1_1_1 and R1_1_2    
    mix$GeneID <- .isoformsToGenes(mix$ID)
    
    # Genes that are defined in the mixture
    geneIDs <- unique(mix$GeneID)

    g <- data.frame(ID=geneIDs,
                    A=rep(0, length(geneIDs)),
                    B=rep(0, length(geneIDs)),
                    Fold=rep(0, length(geneIDs)),
                    LogFold=rep(0, length(geneIDs),
                    Type=rep(NA, length(geneIDs))))

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
        # Type A: Fold change of 1 (negative control)
        # Type B: Fold change other than 1
        #
    

        
    }
    
    # Prefer not to have it as factor variable
    g$ID <- as.character(g$ID)
    
    # Sort by ID so that it'll more easier interpreted
    g <- g[with(g, order(ID)),]
    
    r <- list('isoforms'=data.frame(mix), 'genes'=g)
    class(r) <- c("Mixture")
    r
}
