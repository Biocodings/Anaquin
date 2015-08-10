library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")

files <- c('/Volumes/SSHFS/Sources/QA/scripts/RNA_A1/aligned/accepted_hits.bam',
           '/Volumes/SSHFS/Sources/QA/scripts/RNA_A2/aligned/accepted_hits.bam',
           '/Volumes/SSHFS/Sources/QA/scripts/RNA_A3/aligned/accepted_hits.bam',
           '/Volumes/SSHFS/Sources/QA/scripts/RNA_B1/aligned/accepted_hits.bam',
           '/Volumes/SSHFS/Sources/QA/scripts/RNA_B2/aligned/accepted_hits.bam',
           '/Volumes/SSHFS/Sources/QA/scripts/RNA_B3/aligned/accepted_hits.bam')

bams <- BamFileList(files, yieldSize=2000000)

# Read in the gene model which will be used for counting reads
model <- makeTranscriptDbFromGFF("/Users/tedwong/Sources/QA/data/trans/RNA.v1.gtf", format='gtf')

# Load experimental metadata for the samples
meta <- read.csv(file.path('', "/Users/tedwong/Sources/QA/r/data/experiment.csv"), row.names=1)

# Produces a GRangesList of all the exons grouped by gene
genes <- exonsBy(model, by="gene")

# Fix the duplicat keys
names(bams) <- c("A1.bam", "A2.bam", "A3.bam", "B1.bam", "B2.bam", "B3.bam")

se <- summarizeOverlaps(features=genes, reads=bams, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
	
# The colData slot, so far empty, should contain all the metadata.
colData(se) <- DataFrame(meta)

# We can investigate the resulting SummarizedExperiment by looking at the counts in the assay slot
#head(assay(se))
	


# Build a one-factor model with DESeq2
dds <- DESeqDataSet(se, design = ~condition)

# Run the differential expression
dds <- DESeq(dds)

# Extract the estimated log2 fold changes and p-values for the treated condition
res <- results(dds)

t <- data.frame(id=row.names(res), fold=res$log2FoldChange)





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
