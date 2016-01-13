#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Pool plot is a scatter plot where the variable in the x-axis is categorical. It is a common visualization
# tool for exporing the relationship between sequin groups.
#

plotPool <- function(data, metr,
                     xname = '',
                     yname = '',
                     shouldError = FALSE,
                     cname='Ratio',
                     mix=loadMixture())
{
    require(ggplot2)
    require(RColorBrewer)
    
    stopifnot(!is.null(data))
    stopifnot(!is.null(metr))
    stopifnot(class(data) == 'TransQuin')
    
    seqs <- data$seqs
    
    # Data for the replicates
    reps <- seqs[,c(-1)]
    
    #
    # 1. Calculation the expectation, minimum and maximum
    #
    
    seqs$ave <- rowMeans(reps)
    seqs$max <- apply(reps, 1, max)
    seqs$min <- apply(reps, 1, min)
    
    seqs <- seqs[seqs$ave!=0,]
    
    #
    # 2. Construct the overall abundance
    #
    
    seqs$abund <- NA
    
    for (i in 1:length(1:nrow(seqs)))
    {
        seqs[i,]$abund <- expect(data, id=row.names(seqs[i,]), metr=metr)
    }
    
    # Since the expected abundance varies quite a lot, it's be easier to work on the logarithm scale
    seqs$abund <- log2(seqs$abund)
    
    pal  <- colorRampPalette(brewer.pal(11, "Spectral"))
    cols <- pal(length(unique(seqs$abund)))
    
    y_min <- round(min(log2(seqs$ave)-0.5))
    y_max <- round(max(log2(seqs$ave)+0.5))
    
    #    p <- ggplot(data = seqs, aes(x=log2(X), y=ave, colour=as.factor(abund))) +
    p <- ggplot(data=seqs, aes(x=log2(X), y=log2(ave), colour=abund)) +
        geom_point(size=3)                   +
        xlab(xname)                          +
        ylab(yname)                          +
        labs(colour=cname)                   +
        
        
        #    scale_color_manual(values=cols)      +
        #   scale_fill_manual (values=rev(cols)) +
        scale_colour_gradientn(colours=cols, limits=c(min(seqs$abund), max(seqs$abund)))                  +
        
        scale_x_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        scale_y_continuous(breaks=-5:0, labels=c('1/32','1/16','1/8','1/4','1/2','1'), limits=c(-5,0)) +
        
        #        scale_x_discrete(breaks=c(-2.0,-1.5,-1.0,0), limits=c(-5,0)) +
        #ylim(min(seqs$y) - 0.10, max(seqs$y) + 0.10) +
        #geom_smooth(method = 'lm', formula = y ~ x)  +
        theme_bw()
    
    if (shouldError)
    {
        p <- p + geom_errorbar(aes(ymax=log2(max), ymin=log2(min)), size=0.3, alpha=0.7)
    }
    
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

data <- data[row.names(data)!='R2_63',]
data <- transQuin(seqs=row.names(data), X=data$X, A1=data$A1, A2=data$A2, A3=data$A3, mix=mix)

plotPool(data, metr='gene', shouldError=TRUE, cname='Ratio', xname='Log2 of expected minor/major', yname='Log2 measured minor/major', mix=mix)





