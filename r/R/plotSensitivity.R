#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotSigmoid <- function(data, title, xlab, ylab, showLOA, threshold=0.70, showGuide=FALSE)
{
    require(ggplot2)
    require(reshape2)

    data <- data$seqs

    stopifnot(!is.null(data$input))
    stopifnot(!is.null(data$sensitivity))    
    
    data$f   <- NA
    data$x   <- data$input
    data$y   <- data$sensitivity
    data$grp <- as.factor(round(abs(data$x)))
    
    stopifnot(length(data$x) > 0)
    stopifnot(length(data$x) == length((data$y)) || length(data$x) == nrow((data$y)))
    
    result = tryCatch(
    {
        sigmoid = function(params, x) {
            params[1] / (1 + exp(-params[2] * (x - params[3])))
        }
        
        #
        # Fit a sigmoid curve to the data, which is equivalent to logistic regression.
        #
        
        perf <- min(data[data$y >= 1.00,]$input)
        
        t <- data
        t <- t[t$x <= perf | t$y > 0,]
        x <- t$x
        y <- t$y
        
        fitmodel <- nls(y~a/(1 + exp(-b * (x-c))), start=list(a=1,b=1,c=0))
        params=coef(fitmodel)
        
        data$f <- sigmoid(params, data$x)
    }, error = function(e) {
        showLOA <<- FALSE
    })

    p <- ggplot(data=data, aes(x)) +
                        xlab(xlab) +
                        ylab(ylab) +
                    ggtitle(title) +
                        theme_bw()
    p <- p + geom_point(aes(y=y, colour=grp), size=2.0, alpha=0.5)

    if (!all(is.na(data$f)))
    {
        p <- p + geom_line(aes(y=f, colour="line"), alpha=1.0)
    }
    
    p <- p + theme(axis.title.x=element_text(face='bold', size=12))
    p <- p + theme(axis.title.y=element_text(face='bold', size=12))

    # Limit of assembly
    LOA <- NULL
   
    if (showLOA)
    {
        t <- data[data$f >= threshold,]
        
        if (nrow(t) == 0)
        {
            warning(paste(c('Failed to estimate LOA. The maximum sensitivity is: ', max(data$f))))
        }
        else
        {
            LOA <- round(min(t$input),2)

            label <- 2^LOA
            label <- paste('LOA:', signif(label, 3))
            label <- paste(label, 'attomol/ul')
            
            p <- p + geom_vline(xintercept=LOA, linetype='33')
            p <- p + geom_label(aes(x=min(data$x), y=max(data$y)-0.1), label=label, colour='black', vjust='top', hjust='left', show.legend=FALSE)
        }
    }

    if (!showGuide)
    {
        p <- p + guides(colour=FALSE)
    }
    
    p <- .transformPlot(p)    
    print(p)
    
    return (list(LOA=LOA, fitted=data$f))
}

plotSensitivity <- function(data, title=NULL, xlab=NULL, ylab=NULL, showLOA=TRUE, threshold=0.70)
{
    stopifnot(class(data) == 'Anaquin')
    return (.plotSigmoid(data,
                         title=title,
                         xlab=xlab,
                         ylab=ylab,
                         showLOA=showLOA,
                         threshold=threshold))
}