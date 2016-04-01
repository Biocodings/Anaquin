#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Scatter plot is the most common data visualization tool in Anaquin. It plots the expected concentration
# against the measured concentration. A simple linear regression is fitted on the data.
#

plotScatter <- function(data,
                        alpha = 1.0,
                        shouldLog2 = TRUE,
                        showLegend = FALSE,
                        showStats  = "right",
                        title = '',
                        xname = 'Expected log2 fold change of mixture A and B',
                        yname = 'Measured log2 fold change of mixture A and B')
{
    require(ggplot2)

    data       <- data$seqs
    data$x     <- data$expect
    data$y     <- data$measured
    data$logFC <- abs(log2(data$expect))

    stopifnot(length(data$x) == length((data$y)))
    
    if (length(data$x) == 0)
    {
        warning('Failed to create scatter plot. No data found.')
        return (NULL)
    }
    else if (shouldLog2)
    {
        data$x <- log2(data$x)
        data$y <- log2(data$y)
    }

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
                geom_smooth(method='lm', formula=y ~ x)            +
                labs(colour='Ratio')                               +
                theme_bw()

    x_off <- ifelse(max(data$x) - min(data$x) <= 10, 1.5, 2.0)
    y_off <- ifelse(max(data$y) - min(data$y) <= 10, 0.5, 1.0)

    if (showStats == 'right')
    {
        p <- p + annotate("text", label=lm_eqn(data), x=min(data$x)+x_off, y=max(data$y)-y_off, size=4, colour='black', parse=TRUE)
    }
    else
    {
        p <- p + annotate("text", label=lm_eqn(data), x=min(data$x)+x_off, y=max(data$y)-y_off, size=4, colour='black', parse=TRUE)
    }
    
    p <- p +  theme(axis.title.x=element_text(face='bold', size=15))
    p <- p +  theme(axis.title.y=element_text(face='bold', size=15))

    if (!showLegend)
    {
        p <- p + guides(colour=FALSE)        
    }

    print(p)

    return (list('xname' = xname, 'yname' = yname))
}