#
# This script demostrates how a bedgraph of coverage can be visualized in R.
#

library('Sushi')

# Generated coverage
chrT <- read.csv('/Users/tedwong/Desktop/test.bedgraph', header=FALSE, sep='\t')

# Sequin seqserences
r_seqs <- read.csv('/Users/tedwong/Sources/QA/data/var/AVA017.v032.bed', sep='\t', header=FALSE)

r_seqs <- r_seqs[grep('_R', r_seqs$V4),]

# Reference in the chr21
r_hg38 <- read.csv('/Users/tedwong/Sources/QA/data/var/chr21_random_intervals', sep='\t', header=FALSE)

# Hg38 coverage (chr21)
#hg38 <- read.csv('/Users/tedwong/Desktop/sorted_hg38.bedtools', sep='\t', header=FALSE)
hg38 <- read.csv('/Users/tedwong/Desktop/ira.test', sep='\t', header=FALSE)

colnames(hg38)    <- c('chrom', 'start_', 'end_', 'name', 'pos', 'value')
colnames(chrT)    <- c('chrom', 'start',  'end',  'value')
colnames(r_seqs)  <- c('chrom', 'start',  'end',  'sequin')
colnames(r_hg38)  <- c('chrom', 'start',  'end',  'name')

#
# Convert the non-standard format to bedgraph.
#

hg38$start <- hg38$start_ + (hg38$pos-1)
hg38$end   <- hg38$start_ + (hg38$pos)
    
# Convert to the bedgraph format
hg38 <- subset(hg38, select = c('name', 'start', 'end', 'value'))

# Do we need this?
hg38$name <- 'chr21'

# Draw for each region in the chr21
for (i in 1:nrow(r_hg38))
{
    name  <- as.character(r_hg38[i,]$name)
    start <- as.numeric(r_hg38[i,]$start)
    end   <- as.numeric(r_hg38[i,]$end)
    
    png(filename=paste(paste("/Users/tedwong/Desktop/Coverage/", name, sep=''), '.png', sep=''))
    
    # Construct a density plot for the sequin
    plotBedgraph(hg38, 'chr21', start, end)
    
    labelgenome('chrT', start, end)
    mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
    dev.off()
}

r <- data.frame(id=r_seqs$sequin, coverage=0)

# Draw for each sequin...
for (i in 1:nrow(r_seqs))
{
    id    <- as.character(r_seqs[i,]$sequin)
    start <- as.numeric(r_seqs[i,]$start)
    end   <- as.numeric(r_seqs[i,]$end)

    png(filename=paste(paste("/Users/tedwong/Desktop/Coverage/", id, sep=''), '.png', sep=''))
    
    filter <- chrT[chrT$start >= start & chrT$end <= end,]
    
    # Calculate the average coverage
    r[r$id == id,]$coverage <- mean(filter$value)
    
    # Construct a density plot for the sequin
    plotBedgraph(chrT, 'chrT', start, end)

    labelgenome('chrT', start, end)
    mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
    dev.off()
}

m <- median(r$coverage)
r$fold <- r$coverage/m
r$logFold <- log2(r$fold)
r$diff <- r$coverage - m





