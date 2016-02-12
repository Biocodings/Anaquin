#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotVCoverage <- function(data,
                          alpha = 1.0,
                          title = 'Allele Frequency',
                          xname = 'Expected Allele Fraction (Log2)',
                          yname = 'Coverage (Log2)')
{
    require(ggplot2)

    stopifnot(class(data) == 'VarQuin')

    data  <- data$seqs
    rRead <- data.frame(x=log2(data$expected), y=log2(data$rRead), ratio=data$expected, type=data$type, label='Reference Allele')
    vRead <- data.frame(x=log2(data$expected), y=log2(data$vRead), ratio=data$expected, type=data$type, label='Variant Allele')
    data  <- rbind(rRead, vRead)

    data$ratio <- as.factor(data$ratio)

    # Convert a linear model to string
    lm_eqn <- function(d)
    {
        m <- lm(y ~ x, d);
        eq <- substitute(italic(y) == a + b * italic(x)*','~~italic(r)^2~'='~r2, 
                         list(a  = format(coef(m)[1], digits = 2), 
                              b  = format(coef(m)[2], digits = 2), 
                              r2 = format(summary(m)$r.squared, digits = 3)))
        as.character(as.expression(eq));
    }

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