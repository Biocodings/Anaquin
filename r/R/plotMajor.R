#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotTMajor <- function(data,
                       xname = 'Expected minor fraction (log2)',
                       yname = 'Measured minor fraction (log2)')
{
    require(ggplot2)

    stopifnot(class(data) == 'TransQuin')

    samples <- data$seqs$measured
    row.names(samples) <- row.names(data$seqs)

    # Arithmetic average of all the samples
    samples$avg <- rowMeans(samples)

    # Mapping from isoforms to genes
    samples$gID <- isoformsToGenes(row.names(samples))

    # We can't plot anything that is undefined
    samples <- samples[!is.na(samples$gID),]

    measured <- data.frame(prop=rep(NA, length(unique(samples$gID))), row.names=unique(samples$gID))

    # Expected splicing (minor/major)
    data <- expectForGenes(data, row.names(measured))
    
    for (gID in row.names(measured))
    {
        isos <- samples[samples$gID == gID,]
        stopifnot(nrow(isos) >= 1)
        
        minor <- isos[which.min(isos$avg),]$avg
        major <- isos[which.max(isos$avg),]$avg

        measured[row.names(measured) == gID,] <- minor / major
    }

    data$frac <- measured$prop

    data <- data[data$frac != 0,]
    data <- data[!is.na(data$frac),]

    # Trying to work out the lower-expressed and higher-expressed genes/isoforms.
    data$label <- ifelse(data$expected <= 20, 'Low', 'High')

    # It's be easier to work on the logarithm scale
    data$expected <- log2(data$expected)

    p <- ggplot(data=data, aes(x=log2(efrac), y=log2(frac), color=label)) +
                                          geom_point(size=2) +
                                                 xlab(xname) +
                                                 ylab(yname) +
                                   ggtitle('Minor Fraction') +
                                    labs(colour='Expressed') +        
        scale_x_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        scale_y_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        theme_bw()

    print(p)
}