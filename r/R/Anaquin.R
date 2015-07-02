#  Copyright (C) 2015 - Garvan Institute (Dr Timothey Mercer, Dr Wendy Chen, Ted Wong)
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

LoadMixtures <- function(file)
{
    if (hasArg(file))
    {
        m <- read.csv(file)
    }
    else
    {
        # Nothing specified, load the latest mixture available online
        #mixture <- read.csv(url("http://smallchess.com/Temp/RNA.v1.mix.csv"))
        m <- read.csv('/Users/tedwong/Sources/QA/data/rna/RNA.v4.1.mix', sep='\t')
    }
    
    #
    # The mixture file only gives us the isoforms. We'll need to combine the information for genes
    #
    
    m$GeneID <- as.character(m$ID)
    m$GeneID <- substr(as.character(m$GeneID), 1, nchar(m$GeneID)-2)
    geneIDs  <- unique(m$GeneID)

    # Frame to store data for genes
    g <- data.frame(ID=geneIDs, MixA=rep(0, length(geneIDs)), MixB=rep(0, length(geneIDs)), Fold=rep(0, length(geneIDs)))

    for (gene in geneIDs)
    {
        # Filter a list of sequins for the gene
        sequins <- m[m$GeneID == gene,]

        # Concentration for mixture A at the gene-level
        mixA <- sum(sequins$MixA / sequins$Length)
        
        # Concentration for mixture B at the gene-level
        mixB <- sum(sequins$MixB / sequins$Length)
        
        # Expected fold-change
        fold <- mixB / mixA
        
        g[g$ID == gene,]$MixA <- mixA
        g[g$ID == gene,]$MixB <- mixB        
        g[g$ID == gene,]$Fold <- fold
    }

    g$ID <- as.character(g$ID)
    
    r <- list('data'=data.frame(m), 'genes'=g)
    class(r) <- c("Mixture")
    r
}

print.Mixture <- function(x)
{
    print(head(x$data))
    print(head(x$genes))
}

DESeq2_Analyze <- function(r, m)
{
    # List of genes in the experiment (samples + sequins)
    measured <- rownames(r)

    # List of known genes
    known <- m$genes$ID

    # Filter out the overlapping, the index works only for known
    i <- measured %in% known
  
    # Filter out the overlapping, the index works only for measured
    j <- known %in% measured

    known    <- data.frame(ID=c(known[i]))
    measured <- data.frame(ID=c(measured[j]))
    
    # Known concentration for each sequin detected in the experiment
    known$logFold <- log(m$genes$Fold[i])
    
    # Measured RPKM for each detected sequin
    measured$logFold <- r$log2FoldChange[j]

    stopifnot(nrow(known) == nrow(measured))
    stopifnot(length(known) == length(measured))

    known <- known[with(known, order(ID)),]    
    measured <- measured[with(measured, order(ID)),]
    
    measured$logFold[is.na(measured$logFold)] <- 0
    
    # Fit a simple-linear regression model
    m <- lm(measured$logFold ~ known$logFold)

    # Pearson's correlation
    r <- cor(known$logFold, measured$logFold)

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

EdgeR_Analyze <- function(r, m)
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
    
    # We'll in trouble if they don't match...
    stopifnot(length(genes) == length(sequins))
    
    # Measured RPKM for each detected sequin
    measured <- r$table$logFC[i]
    
    # Concentration for mixture A
    mixA <- m$data$Mix.A[j]
    
    # Concentration for mixture B
    mixB <- m$data$Mix.B[j]
    
    # Known concentration for each detected sequin
    known <- log(mixB / mixA)
    
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
        DESeq2_Analyze(r, m)
    }
    else if (class(r)[1] == 'DGEExact')
    {
        EdgeR_Analyze(r, m)        
    }
    else
    {
        stop('Unknown input. The input must be a result object from DESeq2 or EdgeR')
    }    
}






# 
# 
# 
# #
# # Required: aligned SAM/BAM files
# #
# 
# # Given a list of BAM files, construct an experimental object for DESEq2 and EdgeR
# AnaquinExperiment<- function(files, gtf='/Users/tedwong/Sources/QA/data/rna/RNA.ref.gtf')
# {
# 	bams <- BamFileList(paste('', files, sep='/'), yieldSize=2000000)	
# 
# 	# Read in the gene model which will be used for counting reads
# 	model <- makeTranscriptDbFromGFF(gtf, format='gtf')
# 
# 	# Load experimental metadata for the samples
# 	exp <- read.csv(file.path('', "/Users/tedwong/Sources/QA/r/data/experiment.csv"), row.names=1)
# 	
# 	# Produces a GRangesList of all the exons grouped by gene
# 	genes <- exonsBy(model, by="gene")
# 	
# 	se <- summarizeOverlaps(features=genes, reads=bams, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
# 	
# 	# The colData slot, so far empty, should contain all the metadata.
# 	colData(se) <- DataFrame(exp)
# 
#   se
# 	# We can investigate the resulting SummarizedExperiment by looking at the counts in the assay slot
# 	#head(assay(se))
# 	
# 	# Build a one-factor model with DESeq2
# 	#dds <- DESeqDataSet(se, design = ~condition)
# 	
# 	# Run the differential expression
# 	#dds <- DESeq(dds)
# 	
# 	# Extract the estimated log2 fold changes and p-values for the treated condition
# 	#res <- results(dds)
# 	
# 	# p-values for a particular sequin, one'd expect it be signfiicant due to how the simulation was done
# 	#res['R_1_1',]
# }
# 
# 
# 
# 
# 
# 
# #library("airway")
# #data("airway")
# #se <- airway
# #library("DESeq2")
# #ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
# #d <- DESeq(ddsSE)
# #results(d)
