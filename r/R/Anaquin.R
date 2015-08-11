#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Anaquin is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  Anaquin is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Anaquin If not, see <http://www.gnu.org/licenses/>.

library("GenomicFeatures")

IsoformsToGenes <- function(trans)
{
    trans <- as.character(trans)
    genes <- substr(as.character(trans), 1, nchar(trans)-2)
    genes
}

LoadMixtures <- function()
{
    mix <- read.csv(url('http://anaquin.org/downloads/RNA_4_1.csv'))
    ref <- read.csv(url('http://anaquin.org/downloads/RNA_1.gtf'))
    ref <- makeTranscriptDbFromGFF(file = '/Users/tedwong/Sources/QA/data/trans/RNA.v1.gtf', format = "gtf")
    ref <- unique(IsoformsToGenes(transcripts(ref)$tx_name))
    
    # Eg: R1_1 for R1_1_1 and R1_1_2    
    mix$geneID <- IsoformsToGenes(mix$ID)
    
    # Genes that are defined in the mixture
    geneIDs <- unique(mix$geneID)
    
    # Genes that are defined in reference but not reference and therfore must be ignored as there's no concentration
    ignored <- ref[!(ref %in% geneIDs)]
    
    combined <- c(geneIDs, ignored)
    
    g <- data.frame(ID=combined,
                    a=rep(0, length(combined)),
                    b=rep(0, length(combined)),
                    fold=rep(0, length(combined)),
                    logFold = rep(0, length(combined)))
    
    for (id in geneIDs)
    {
        seqs <- mix[mix$geneID == id,]
        
        #
        # Calculate the expected abundance
        #
        
        g[g$ID == id,]$a <- sum(seqs$MixA)
        g[g$ID == id,]$b <- sum(seqs$MixB)
        
        #
        # Calculate the expected fold-ratio
        #
        
        g[g$ID == id,]$fold    <- g[g$ID == id,]$b / g[g$ID == id,]$a
        g[g$ID == id,]$logFold <- log2(g[g$ID == id,]$fold)
    }
    
    #
    # We simply can't assume that sequins defined in the mixture and reference are identical.
    # Here, we find out those sequins that are not defined in the mixture.
    #
    
    for (id in ignored)
    {
        g[g$ID == id,]$a       <- 'NA'
        g[g$ID == id,]$b       <- 'NA'
        g[g$ID == id,]$fold    <- 'NA'
        g[g$ID == id,]$logFold <- 'NA'
    }
    
    # Prefer not to have it as factor variable
    g$ID <- as.character(g$ID)
    
    # Sort by ID so that it'll more easier interpreted
    g <- g[with(g, order(ID)),]
    
    r <- list('data'=data.frame(mix), 'genes'=g)
    class(r) <- c("Mixture")
    r
}

print.Mixture <- function(x)
{
    print(head(x$data))
    print(head(x$genes))
}

DESeq2 <- function(r, mix)
{
    # List of known genes
    known <- as.character(mix$genes$ID)
    
    # Sequins that are detected in the experiment    
    detected <- rownames(r) %in% known
    
    # Filter out to only the rows with sequins
    r <- r[detected,]
    
    # Create a data-frame for each sequin defined in the mixture and reference, whether it's been detected
    d <- data.frame(id=known, known=rep(NaN, length(known)), measured=rep(NaN, length(known)))
    
    # For each sequin detected in the experiment
    for (id in rownames(r))
    {
        d[d$id==id,]$known    <- mix$genes[mix$genes$ID==id,]$logFold
        d[d$id==id,]$measured <- r[id,]$log2FoldChange
    }
    
    #
    # Fit a linear model for sequins that are detected in the experiment.
    #
    
    d_ <- d[is.finite(d$measured),]
    
    # Fit a simple-linear regression model
    m <- lm(d_$known ~ d_$measured)
    
    # Pearson's correlation
    r <- cor(as.numeric(d_$known), as.numeric(d_$measured))
    
    # Coefficients of determination
    r2 <- summary(m)$r.squared
    
    # Regression slope
    slope <- coef(m)["known"]
    
    # Generate a linear plot of the relationship
    plot(d_$known, d_$measured)
}

EdgeR <- function(r, m)
{
    # List of genes for the experiment (sequins + samples)
    genes <- rownames(r)
    
    # List of all sequins
    sequins <- c(as.vector(m$data$ID))
    
    # Filter out to only rows for sequins, the index works only on genes
    i <- genes %in% sequins
    
    # Filter out undetected sequins, the index works only on sequins
    j <- sequins %in% genes
    
    genes   <- genes[i]
    sequins <- sequins[j]
    
    stopifnot(length(genes) == length(sequins))
    
    # Measured RPKM for each detected sequin
    measured <- r$table$logFC[i]
    
    # Concentration for mixture A
    mixA <- m$data$Mix.A[j]
    
    # Concentration for mixture B
    mixB <- m$data$Mix.B[j]
    
    # Known concentration for each detected sequin
    known <- log2(mixB / mixA)
    
    # We'll again in trouble if they don't match...
    stopifnot(length(known) == length(measured))
    
    # Fit a simple-linear regression model
    m <- lm(measured~known)
    
    # Pearson's correlation
    r <- cor(known, measured)
    
    # Coefficients of determination
    r2 <- summary(m)$r.squared
    
    # Regression slope
    slope <- coef(m)["known"]
    
    # Generate a linear plot of the relationship
    plot(known, measured)
    
    r <- list(data=data.frame(known, measured), r=r, r2=r2, slope=slope)
    class(r) <- c("Anaquin")
    r 
}

Anaquin <- function(r, m=LoadMixtures())
{
    if (class(r) == 'DESeqResults')
    {
        DESeq2(r, m)
    }
    else if (class(r)[1] == 'DGEExact')
    {
        EdgeR(r, m)        
    }
    else
    {
        stop('Unknown input. The input must be a result object from DESeq2 or EdgeR')
    }    
}
