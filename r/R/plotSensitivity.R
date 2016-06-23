#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotSigmoid <- function(data, threshold, title='', xlab='', ylab='', showLOQ=TRUE, showGuide=FALSE)
{
    require(ggplot2)
    require(reshape2)

    data <- data$seqs

    stopifnot(!is.null(data$expected))
    stopifnot(!is.null(data$measured))    
    
    data$f   <- NA
    data$x   <- data$expected
    data$y   <- data$measured
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
            
            perf <- min(data[data$y >= 1.00,]$expected)
            
            t <- data
            t <- t[t$x <= perf | t$y > 0,]
            x <- t$x
            y <- t$y
            
            fitmodel <- nls(y~a/(1 + exp(-b * (x-c))), start=list(a=1,b=1,c=0))
            params=coef(fitmodel)
            data$f <- sigmoid(params,data$x)
        }, error = function(e)
        {
            showLOQ <<- FALSE
            #warning(e)
        })

    p <- ggplot(data=data, aes(x)) +
                        xlab(xlab) +
                        ylab(ylab) +
                    ggtitle(title) +
                        theme_bw()
    p <- p + geom_point(aes(y=y, colour=grp), alpha=1.0)
    
    if (!is.na(data$f))
    {
        p <- p + geom_line(aes(y=f, colour="line"), alpha=1.0)
    }
    
    p <- p + theme(axis.title.x=element_text(face='bold', size=12))
    p <- p + theme(axis.title.y=element_text(face='bold', size=12))
    
    if (showLOQ)
    {
        r <- round(min(data[data$y >= threshold,]$expected),2)
        
        # We can assume the break-point is on the log2-scale. Let's convert it back.
        label <- 2^r

        label <- paste('LOA:', signif(label, 3))
        label <- paste(label, 'attomol/ul')

        p <- p + geom_vline(xintercept=r, linetype='33')
        p <- p + geom_label(aes(x=min(data$x), y=max(data$y)), label=label, colour='black', vjust='top', hjust='left', show.legend=FALSE)
    }

    if (!showGuide)
    {
        p <- p + guides(colour=FALSE)
    }
    
    p <- .transformPlot(p)    
    print(p)
}

plotSensitivity <- function(data, title, xlab, ylab, showLOQ=TRUE, threshold=0.98)
{
    stopifnot(class(data) == 'Anaquin')
    .plotSigmoid(data, title=title, xlab=xlab, ylab=ylab, showLOQ=showLOQ, threshold=threshold)    
}