
dir <- '/Users/tedwong/Sources/QA/scripts/test'
files <- c('/Users/tedwong/Sources/QA/scripts/test/mix_a.bam', '/Users/tedwong/Sources/QA/scripts/test/mix_b.bam')

library("Rsamtools")
bamfiles <- BamFileList(files, yieldSize=2000000)

library("GenomicFeatures")

gtffile <- '/Users/tedwong/Sources/QA/scripts/test/standards.gtf'
(txdb <- makeTranscriptDbFromGFF(gtffile, format="gtf"))

(genes <- exonsBy(txdb, by="gene"))

library("GenomicAlignments")
se <- summarizeOverlaps(features=genes, reads=bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )

# Meta-data for the experiment
sampleTable <- read.csv('/Users/tedwong/Sources/QA/scripts/test/experiment.csv', row.names=1)

colData(se) <- DataFrame(sampleTable)

library("DESeq2")
dds <- DESeqDataSet(se, design = ~Condition)

# Need replicates, just put them as new library in the CSV!!!!





