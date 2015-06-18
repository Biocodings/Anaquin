# #  Copyright (C) 2015 - Garvan Institute (Dr Timothey Mercer, Dr Wendy Chen, Ted Wong)
# #
# #  Spike is free software: you can redistribute it and/or modify
# #  it under the terms of the GNU General Public License as published by
# #  the Free Software Foundation, either version 2 of the License, or
# #  (at your option) any later version.
# #
# #  Spike is distributed in the hope that it will be useful,
# #  but WITHOUT ANY WARRANTY; without even the implied warranty of
# #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# #  GNU General Public License for more details.
# #
# #  You should have received a copy of the GNU General Public License
# #  along with Spike. If not, see <http://www.gnu.org/licenses/>.
# 
# #
# # Demostrates common usage of sequin analysis with DESeq2 and EdgeR. We will follow a simulated data-set
# # for an experiment. There is a single factor in the experiment, the second group has ten fold-change
# # relative to the first group.
# #
# 
# #
# # Required: aligned SAM/BAM files
# #
# 
# # Given a list of BAM files, construct an experimental object for DESEq2 and EdgeR
# AnaquinExperiment<- function(gtf='data/rna/standards.gtf')
# {
# 	bams <- BamFileList(paste('/Users/tedwong/Sources/QA/r/data/', list.files(path = '/Users/tedwong/Sources/QA/r/data/', pattern = "\\.bam$"), sep='/'), yieldSize=2000000)
# 	
# 	# Read in the gene model which will be used for counting reads
# 	model <- makeTranscriptDbFromGFF('/Users/tedwong/Sources/QA/data/rna/RNA.ref.gtf', format='gtf')
# 
# 	# Load experimental metadata for the samples
# 	exp <- read.csv(file.path('/Users/tedwong/Sources/QA/r/data/', "experiment.csv"), row.names=1)
# 	
# 	# Produces a GRangesList of all the exons grouped by gene
# 	genes <- exonsBy(model, by="gene")
# 	
# 	se <- summarizeOverlaps(features=genes, reads=bams, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)	
# 
# 	# The colData slot, so far empty, should contain all the metadata.
# 	colData(se) <- DataFrame(exp)
#     
#     se
# }
# 
# library("DESeq2")
# library("Rsamtools")
# library("GenomicFeatures")
# library("GenomicAlignments")
# 
# 
# 
# # We can investigate the resulting SummarizedExperiment by looking at the counts in the assay slot
# #head(assay(se))
# 
# #
# # R_1_4 is missing in the experiment as expected because it was skipped in the simulation.
# #
# 
# # The colData slot, so far empty, should contain all the metadata.
# #colData(se) <- DataFrame(exp)
# 
# # Build a one-factor model with DESeq2
# dds <- DESeqDataSet(se, design = ~condition)
# 
# # Run the differential expression
# dds <- DESeq(dds)
# 
# # Extract the estimated log2 fold changes and p-values for the treated condition
# res <- results(dds)
# 
# # p-values for a particular sequin, one'd expect it be signfiicant due to how the simulation was done
# res['R_1_1',]
# 
# # What does Anaquin do? Anaquin does differential expression here and compare the log-fold changes
# #
# # The column log2FoldChange!!!!
# 
# 
# 
# 
# 
# 
# 
# 
# 
