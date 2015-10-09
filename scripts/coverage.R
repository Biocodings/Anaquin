

x <- read.csv('/home/tedwon/Coverage/chrT.bedgraph', header=FALSE, sep='\t')
y <- read.csv('/home/tedwon/Coverage/hg38.bedgraph', header=FALSE, sep='\t')



x <- read.csv('/Users/tedwong/Sources/chrT.bedgraph', header=FALSE, sep='\t')

chrT <- data.frame(coverage=x$V4)
hg38 <- data.frame(coverage=y$V4)
hg38$coverage = hg38$coverage*2

chrT$veg <- 'chrT'
hg38$veg <- 'hg38'

vegLengths <- rbind(chrT, hg38)
vegLengths <- rbind(chrT)

ggplot(vegLengths, aes(coverage, fill = veg)) + geom_density(alpha = 0.2)
ggsave('TimMercer.pdf')







#p1 <- hist(chrT$data)
p2 <- hist(hg38$data)

#plot(p1, col=rgb(0,0,1,1/4))
#plot(p2, col=rgb(1,0,1,1/4), add=T)
plot(p2, col=rgb(1,0,1,1/4))

dev.off()

# Try to do the same to hg38 and show the graphs together...



    # run('grep -v chrT aligned.sam > hg38.sam')
    #   run('grep chrT aligned.sam > chrT.sam')
    #  run('samtools view -bS chrT.sam > chrT.bam')
    
    # run('wc hg38.sam') # How many reads? 150586265
    #run('wc chrT.sam') # How many reads? 150586265
    
    #
    # It's probably too time-consuming to construct a distribution for the genome. Why not just take a sufficent random samples?
    #
    
    
    
    #run('samtools view -s 0.1727372 -b chrT.bam > sampled.bam')

    #    run('samtools sort chrT.bam aligned_chrT')
    #   run('samtools sort sampled.bam aligned_sampled')
    #  run('bedtools genomecov -bg -ibam aligned_chrT.bam > aligned_chrT.bedgraph')
    # run('bedtools genomecov -bg -ibam aligned_sampled.bam > aligned_sampled.bedgraph')

    #
    # Draw each individual bedgraph on the same plot...
    #

    #    run('grep chrT aligned.bedgraph > chrT.bedgraph')
    #    run('grep -v chrT aligned.bedgraph > hg38.bedgraph')