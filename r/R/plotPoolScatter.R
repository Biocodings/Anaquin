#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Scatter plot is the most common data visualization tool in Anaquin. It plots the expected concentration
# against the measured concentration. A simple linear regression is fitted on the data.
#

plotPoolScatter <- function(data,
                            alpha = 1.0,
                            shouldLog2 = TRUE,
                            shouldHideLegend = TRUE,
                            title = '',
                            xname = 'Expected log2 fold change of mixture A and B',
                            yname = 'Measured log2 fold change of mixture A and B')
{
    require(ggplot2)

    data       <- data$seqs
    data$x     <- data$expected
    data$A1    <- data$A1
    data$A2    <- data$A2
    data$A3    <- data$A3    
    data$logFC <- abs(log2(data$expected))
    
    if (shouldLog2)
    {
        data$x <- log2(data$x)
        data$A1 <- log2(data$A1)
        data$A2 <- log2(data$A2)
        data$A3 <- log2(data$A3)
    }

    data$y    <- rowMeans(data[c(2,3,4)])
    data$sd   <- apply(data[c(2,3,4)], 1, sd)
    data$ymax <- data$y + data$sd
    data$ymin <- data$y - data$sd
    
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

    data$logFC <- as.factor(data$logFC)
    
    p <- ggplot(data=data, aes(x=x, y=y)) +
                              xlab(xname) +
                              ylab(yname) +
                           ggtitle(title) +
                geom_point(aes(colour=logFC), size=2, alpha=alpha) + 
                #xlim(min(data$x) - 0.10, max(data$x) + 0.10)       +
                #ylim(min(data$y) - 0.10, max(data$y) + 0.10)       +
                geom_smooth(method='lm', formula=y ~ x)            +
                annotate("text", label=lm_eqn(data), x=0, y=max(data$y)-1, size=5, colour='black', parse=TRUE) +
                theme_bw()

    p <- p +  theme(axis.title.x=element_text(face='bold', size=15))
    p <- p +  theme(axis.title.y=element_text(face='bold', size=15))
    p <- p + geom_errorbar(aes(ymax=ymax, ymin=ymin), size=0.3, alpha=0.7)

    if (shouldHideLegend)
    {
        p <- p + guides(colour=FALSE)        
    }
    
    print(p)

    return (list('xname' = xname, 'yname' = yname))
}