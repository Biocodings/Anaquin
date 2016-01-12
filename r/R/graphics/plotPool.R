#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Pool plot is a scatter plot where the variable in the x-axis is categorical. It is a common visualization
# tool for exporing the relationship sequin groups. To improve the visibility, the coordinatives are not
# drawn by scale.
#

plotPool <- function(data,
                     cname = 'LogFC', 
                     xname = 'Expected log2 fold change of minor/major',
                     yname = 'Measured log2 fold change of minor/major')
{
    require(ggplot2)
    require(RColorBrewer)
    
    data <- data[!is.na(data$A1),]
    data <- data[!is.na(data$A2),]
    data <- data[!is.na(data$A3),]

    # Data for the replicates
    reps <- data[,c(-1)]
    
    data$ave <- rowMeans(reps)
    data$max <- apply(reps, 1, max)
    data$min <- apply(reps, 1, min)
    
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
    
    pal <- colorRampPalette(rev(brewer.pal(length(unique(data$X)), "Spectral")))

    data$x <- data$X
    data$y <- data$ave
    data$ratio <- as.factor(log2(data$x))
    
    p <- ggplot(data = data, aes(x = log2(x), y = y, colour=ratio))    +
                    xlab(xname)                                  +
        labs(colour=cname) +
                    ylab(yname)                                  +
                    geom_point(size = 3)                         +
        
        
  #      geom_errorbar(aes(ymax = max, ymin = min,
   #                       colour = ratio), size = 0.3, alpha = 0.4)              +
        
#        scale_x_discrete(breaks=c(-2.0,-1.5,-1.0,0.5,0), limits=c(-1.5:0)) +
                    #xlim(min(data$x) - 0.10, max(data$x) + 0.10) +
                    #ylim(min(data$y) - 0.10, max(data$y) + 0.10) +
                    #geom_smooth(method = 'lm', formula = y ~ x)  +
                    #annotate("text", label = lm_eqn(data), x = 0, y = max(data$y), size = 6, colour = 'black', parse=TRUE) +
                    #scale_colour_gradientn(colours = pal(20), limits=c(min(data$logFC), max(data$logFC)))                  +
                    theme_bw()
    print(p)
    
    return (list('xname' = xname, 'yname' = yname))
}




data <- read.csv('/Users/tedwong/Desktop/K_RMXA_minIsoform_full_3reps_TED.csv', row.names=1)
data <- data[-c(1),]
colnames(data) <- c('X', 'A1', 'A2', 'A3')

data <- data[data$A1!='#DIV/0!',]
data <- data[data$A2!='#DIV/0!',]
data <- data[data$A3!='#DIV/0!',]

data$X  <- as.numeric(as.character(data$X))
data$A1 <- as.numeric(as.character(data$A1))
data$A2 <- as.numeric(as.character(data$A2))
data$A3 <- as.numeric(as.character(data$A3))


plotPool(data)
