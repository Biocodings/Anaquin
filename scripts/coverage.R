
x <- read.csv('/home/tedwon/Coverage/chrT.bedgraph', header=FALSE, sep='\t')
y <- read.csv('/home/tedwon/Coverage/hg38.bedgraph', header=FALSE, sep='\t')

chrT <- data.frame(coverage=x$V4)
hg38 <- data.frame(coverage=y$V4)
hg38$coverage = hg38$coverage*2

chrT$veg <- 'chrT'
hg38$veg <- 'hg38'

#vegLengths <- rbind(chrT, hg38)
vegLengths <- rbind(hg38)

ggplot(vegLengths, aes(coverage, fill = veg)) + geom_density(alpha = 0.2)
ggsave('TimMercer.pdf')







#p1 <- hist(chrT$data)
p2 <- hist(hg38$data)

#plot(p1, col=rgb(0,0,1,1/4))
#plot(p2, col=rgb(1,0,1,1/4), add=T)
plot(p2, col=rgb(1,0,1,1/4))

dev.off()

# Try to do the same to hg38 and show the graphs together...

