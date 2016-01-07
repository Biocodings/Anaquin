#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Scatter plot is the most common data visualization tool in Anaquin. It plots the expected concentration
# defined by a mixture with the measurement.
#

plotScatter <- function(data,
                        shouldLog2 = TRUE,
                        xname = 'Expected log2 fold change of mixture A and B',
                        yname = 'Measured log2 fold change of mixture A and B')
{
    require(ggplot2)

    stopifnot(class(data) == 'ExpectMeasured')

    data   <- data[['data']]
    data$x <- data$expected
    data$y <- data$measured
    
    if (shouldLog2)
    {
        data$x <- log2(data$x)
        data$y <- log2(data$y)
    }

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
    
    p <- ggplot(data = data, aes(x = x, y = y))                 +
                    theme_bw()                                  +
                    xlab(xname)                                 +
                    ylab(yname)                                 +
                    geom_point()                                +
                    ggtitle('')                                 +
                    xlim(min(data$x) - 1, max(data$x) + 1)      +
                    ylim(min(data$y) - 1, max(data$y) + 1)      +
                    geom_smooth(method = 'lm', formula = y ~ x) +
                    geom_text(x = 0, y = max(data$y) + 1, label = lm_eqn(data), parse = TRUE)
    print(p)
}