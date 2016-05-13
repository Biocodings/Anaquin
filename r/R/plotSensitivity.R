#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotSensitivity <- function(data, limit, title='', xlab='', ylab='', showLimit=TRUE, showGuide=FALSE)
{
    require(ggplot2)
    require(reshape2)

    data     <- data$seqs
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
        showLimit <<- FALSE
        warning(e)
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

    print (showLimit)
    
    if (showLimit)
    {
        r <- min(data[data$y >= limit,]$expected)
        p <- p + geom_vline(xintercept=c(r), linetype="dotted")
        p <- p + geom_label(aes(x=r, y=0.30, label=paste('LOA', r)), colour = "black", show.legend=FALSE)
        #data[data$x > loa & data$y == 0,]
    }
    
    if (!showGuide)
    {
        p <- p + guides(colour=FALSE)
    }
    
    print(p)
}

plotSensitivity.FusQuin <- function(data, title, xlab, ylab, limit)
{
    .plotSensitivity(data, title=title, xlab=xlab, ylab=ylab, limit=limit)
}

plotSensitivity.MetaQuin <- function(data, title, xlab, ylab, limit)
{
    .plotSensitivity(data, title=title, xlab=xlab, ylab=ylab, limit=limit)
}

plotSensitivity.TransQuin <- function(data, title, xlab, ylab, limit)
{
    .plotSensitivity(data, title=title, xlab=xlab, ylab=ylab, limit=limit)
}

plotSensitivity <- function(data, title=NULL, xlab=NULL, ylab=NULL, limit=0.98)
{
    stopifnot(class(data) == 'FusQuin'   |
              class(data) == 'TransQuin' |
              class(data) == 'MetaQuin')

    if (class(data) == 'FusQuin')   { plotSensitivity.FusQuin(data, title=title, xlab=xlab, ylab=ylab, limit=limit) } 
    if (class(data) == 'TransQuin') { plotSensitivity.TransQuin(data, title=title, xlab=xlab, ylab=ylab, limit=limit) } 
    if (class(data) == 'MetaQuin')  { plotSensitivity.MetaQuin(data, title=title, xlab=xlab, ylab=ylab, limit=limit) } 
}