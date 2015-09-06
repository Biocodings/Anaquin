library('Sushi')

d <- read.csv('/Users/tedwong/Desktop/sorted.bedgraph', header=FALSE, sep='\t')
c <- read.csv('/Users/tedwong/Desktop/Untitled.csv', header=FALSE)

for (i in 1:nrow(c))
{
    id <- as.character(c[i,]$V1)
    start <- as.numeric(c[i,]$V2)
    end <- as.numeric(c[i,]$V3)

    png(filename=paste(paste("/Users/tedwong/Desktop/", id, sep=''), '.png', sep=''))
    plotBedgraph(d, 'chrT', start-1000, end+1000)
    labelgenome('chrT', start, end)
    mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
    dev.off()
}



