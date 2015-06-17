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
        mixture <- read.csv(file)
    }
    else
    {
        # Nothing specified, load the latest mixture file available online
        mixture <- read.csv(url("http://smallchess.com/Temp/RNA.v1.mix.csv"))
    }

    version <- unique(mixture[["Version"]])

    r <- c(mixture=mixture, version=version)
    #class(r) <- c("Mixture")
    r
}

#print.Mixture <- function(x)
#{
#    print(paste('Version:', x$version))
#    print(x$mixture)
#}

RNAAnalyze <- function(r)
{
    # Load the default and latest mixtures
    m <- LoadMixtures()
  
    c <- assays(se)$counts

    sequins <- c(as.vector(m$mixture.ID))
  
    i <- row.names(c) %in% sequins
  
    j <- sequins %in% row.names(c)
  
    # Filter out the rows for sequins (other rows are assumed to be the experiment)
    rows <- c[i,]
  
    # Measured RPKM for each detected sequin
    observed <- res$log2FoldChange[i]
  
    # Concentration for mixture A
    mixA <- m$mixture.Mix.A[j]
  
    # Concentration for mixture B
    mixB <- m$mixture.Mix.B[j]
  
    # Known concentration for each detected sequin
    known <- log(mixB / mixA)
}









#
# Required: aligned SAM/BAM files
#

# Given a list of BAM files, construct an experimental object for DESEq2 and EdgeR
AnaquinExperiment<- function(files, gtf='/Users/tedwong/Sources/QA/data/rna/RNA.ref.gtf')
{
	bams <- BamFileList(paste('', files, sep='/'), yieldSize=2000000)	

	# Read in the gene model which will be used for counting reads
	model <- makeTranscriptDbFromGFF(gtf, format='gtf')

	# Load experimental metadata for the samples
	exp <- read.csv(file.path('', "/Users/tedwong/Sources/QA/r/data/experiment.csv"), row.names=1)
	
	# Produces a GRangesList of all the exons grouped by gene
	genes <- exonsBy(model, by="gene")
	
	se <- summarizeOverlaps(features=genes, reads=bams, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
	
	# The colData slot, so far empty, should contain all the metadata.
	colData(se) <- DataFrame(exp)

  se
	# We can investigate the resulting SummarizedExperiment by looking at the counts in the assay slot
	#head(assay(se))
	
	# Build a one-factor model with DESeq2
	#dds <- DESeqDataSet(se, design = ~condition)
	
	# Run the differential expression
	#dds <- DESeq(dds)
	
	# Extract the estimated log2 fold changes and p-values for the treated condition
	#res <- results(dds)
	
	# p-values for a particular sequin, one'd expect it be signfiicant due to how the simulation was done
	#res['R_1_1',]
}






#library("airway")
#data("airway")
#se <- airway
#library("DESeq2")
#ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
#d <- DESeq(ddsSE)
#results(d)
