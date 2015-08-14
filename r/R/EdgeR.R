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

library('edgeR')

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




# Factor level
group <- factor(c(1,1,1,2,2,2))

y <- DGEList(counts=assay(se), group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
r <- exactTest(y)




