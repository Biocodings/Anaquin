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



p <- ggplot(data = r, aes(x = r$known, y = r$measured))
p <- p + xlab('Log2 known coverage')
p <- p + ylab('Log2 measured log-fold')
p <- p + geom_point()
p <- p + theme_classic()
p <- p + geom_smooth(method = "lm")
print(p)


