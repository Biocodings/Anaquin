#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotAllele <- function(data)
{
    .plotExpress(data, title='Allele frequency',
                        xlab='Expected allele frequency (log2)',
                        ylab='Measured allele frequency (log2)',
                     showLOQ=FALSE,
                   showStats='left')
}

plotAlleleReads <- function(data,
                            alpha = 1.0,
                            title = 'Allele Frequency vs Read count',
                             xlab = 'Allele Frequency (Log2)',
                             ylab = 'Observed read count (Log2)')
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