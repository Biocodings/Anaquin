#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotCoverage <- function(data, ylab='Sequence Coverage', seqs=NULL, minLen=NULL, maxLen=NULL)
{
    require(plyr)
    require(Sushi)
    require(rtracklayer)

    #
    # 1. Sequin annotation in BED format
    # 2. Generated coverage in BEDGRAPH format
    #
    
    stopifnot(!is.null(data$annot))
    stopifnot(!is.null(data$bedgr))    

    # Read the reference annotation
    ref <- import.bed(con=data$annot)

    # Read the generated coverage
    data <- read.csv(data$bedgr, header=FALSE, sep='\t')

    colnames(data) <- c('chrom', 'start',  'end',  'value')

    for (i in 1:length(ref))
    {
        r <- ref[i,]
        
        #
        # Do we want to filter to a particular set of sequins?
        #
        
        if (is.null(seqs) || seq_$name %in% seqs)
        {
            startL <- start(r)
            endL   <- end(r)
            
            #
            # Do we want to filter by minimum/maximum length?
            #
            
            if ((is.null(minLen) || startL >= minLen) && (is.null(maxLen) || endL <= maxLen))
            {
                Sushi::plotBedgraph(data, 'chrT', start(r), end(r), transparency=0.50, color='#ADC2E6', xlab=VarQuin.genes(r$name), ylab=ylab)

                ticks  <- 5
                range  <- c(min(data[data$start >= start(r) & data$end <= end(r),]$value),
                            max(data[data$start >= start(r) & data$end <= end(r),]$value))
                
                scaled <- range / ticks
                scaled <- round_any(scaled, 100, f = ceiling)
                
                axis(side=2, at=seq(0, ticks * scaled[2], by=scaled[2]))
                labelgenome('chrIS', start(r), end(r), n=4, scale="Mb")
            }
        }
    }
}