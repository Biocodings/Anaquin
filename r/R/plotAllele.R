#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotAllele <- function(data)
{
    .plotExpress(data, title='Expected allele frequency vs Measured allele frequency',
                       xname='Expected allele frequency (log2)',
                       yname='Measured allele frequency (log2)',
                   showStats='left',
                  showLegend=TRUE)
}

plotAlleleReads <- function(data,
                            alpha = 1.0,
                            title = 'Allele Frequency vs Read count',
                            xname = 'Allele Frequency (Log2)',
                            yname = 'Observed read count (Log2)')
{
    require(ggplot2)

    stopifnot(class(data) == 'VarQuin')

    data  <- data$seqs
    rRead <- data.frame(x=log2(data$expect), y=log2(data$rRead), ratio=data$expect, type=data$type, label='Reference Allele')
    vRead <- data.frame(x=log2(data$expect), y=log2(data$vRead), ratio=data$expect, type=data$type, label='Variant Allele')
    data  <- rbind(rRead, vRead)

    data$ratio <- as.factor(data$ratio)

    p <- ggplot(data=data, aes(x=x, y=y)) +
                              xlab(xname) +
                              ylab(yname) +
                           ggtitle(title) +
                             geom_point() +
                        scale_x_reverse() +
                         labs(colour='')  +
            geom_point(aes(colour=label), size=2, alpha=alpha) +
                theme_bw()

    p <- p +  theme(axis.title.x=element_text(face='bold', size=15))
    p <- p +  theme(axis.title.y=element_text(face='bold', size=15))

    print(p)

    return (list('xname' = xname, 'yname' = yname))
}