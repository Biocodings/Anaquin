
# Read in for the base coverage...
d <- read.csv('/Users/tedwong/Sources/QA/chrT_sorted.bedgraph', header=FALSE, sep='\t')

dd <- d$V4

hist(dd)

# Try to do the same to hg38 and show the graphs together...

