#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotDensity <- function(data, seqIDs=NULL, minBase=NULL, maxBase=NULL)
{
    require(plyr)
    require(Sushi)
    require(rtracklayer)

    # Read the reference annotation
    ref <- import.bed(con=data$annot)

    data <- read.csv(data$bedgr, header=FALSE, sep='\t')
    colnames(data) <- c('chrom', 'start',  'end',  'value')

    for (i in 1:length(ref))
    {
        r <- ref[i,]
        
        if (is.null(seqIDs) || seq_$name %in% seqIDs)
        {
            startL <- start(r)
            endL   <- end(r)
            
            if ((is.null(minBase) || startL >= minBase) && (is.null(maxBase) || endL <= maxBase))
            {
                Sushi::plotBedgraph(data, 'chrT', start(r), end(r), transparency=0.50, color='#ADC2E6', xlab=r$name, ylab='Coverage')

                ticks  <- 5
                range  <- c(min(data[data$start >= start(r) & data$end <= end(r),]$value),
                            max(data[data$start >= start(r) & data$end <= end(r),]$value))
                
                scaled <- range / ticks
                scaled <- round_any(scaled, 100, f = ceiling)
                
                axis(side=2, at=seq(0, ticks * scaled[2], by=scaled[2]))
            }
        }
    }
}