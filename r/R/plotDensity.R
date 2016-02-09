#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotDensity <- function(file, seqIDs = NULL, minBase = NULL, maxBase = NULL)
{
    require(Sushi)
    require(plyr)
    require(rtracklayer)

    data <- read.csv(file, header=FALSE, sep='\t')
    colnames(data) <- c('chrom', 'start',  'end',  'value')

    # Read the reference annotation
    ref <- import.bed(con='/Users/tedwong/Sources/QA/data/VARQuin/AVA017.v032.bed')

    # Loop over all sequins defined...
    for (i in 1:length(ref))
    {
        seq_ <- ref[i,]
        
        if (is.null(seqIDs) || seq_$name %in% seqIDs)
        {
            startL <- start(seq_)
            endL   <- end(seq_)
            
            if ((is.null(minBase) || startL >= minBase) && (is.null(maxBase) || endL <= maxBase))
            {
                Sushi::plotBedgraph(data, 'chrT', start(seq_), end(seq_), transparency=0.50, color='#ADC2E6', xlab=seq_$name, ylab='Coverage')
                #Sushi::plotBedgraph(data, 'chrT', 6535259, 6536017, transparency=0.50, color='#ADC2E6', xlab=seq_$name, ylab='Coverage')
                
                ticks  <- 5
                
                range  <- c(min(data[data$start >= start(seq_) & data$end <= end(seq_),]$value),
                            max(data[data$start >= start(seq_) & data$end <= end(seq_),]$value))
                
                scaled <- range / ticks
                scaled <- round_any(scaled, 100, f = ceiling)
                
                axis(side=2, at=seq(0, ticks * scaled[2], by=scaled[2]))
            }
        }
    }
}