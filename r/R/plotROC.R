#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotROC <- function(data, refRats, ...)
{
    data <- data$seqs
    
    stopifnot(!is.null(refRats))
    stopifnot(!is.null(data$score))
    stopifnot(!is.null(data$label))
    stopifnot(!is.null(data$input))    
    
    x <- list(...)

    if (is.null(x$title))      { x$title      <- NULL    }
    if (is.null(x$legTitle))   { x$legTitle   <- 'Ratio' }
    if (is.null(x$showLegend)) { x$showLegend <- TRUE    }

    # This is the sequin groups
    data$ratio <- abs(round(data$input))

    data <- data[!is.na(data$score),]
    data <- data[data$label=='TP' | data$label=='FP',]
    data <- data[, order(names(data))]
    data <- data.frame(label=data$label, score=data$score, ratio=data$ratio)

    ROCs   <- NULL
    AUCs   <- NULL
    ratios <- sort(data$ratio)

    # Query ratios (not including the references)
    uniqs <- unique(ratios)
    
    # Remove the reference ratio (if any)
    uniqs <- uniqs[uniqs != refRats & !(uniqs %in% refRats)]
    
    stopifnot(length(uniqs) > 0)

    #
    # The reference ratios should have the same length as the queries. If there is only
    # a single reference, maybe we can replicate it to the query length?
    #
    
    if (length(refRats) == 1 && length(refRats) != length(uniqs))
    {
        refRats <- rep(refRats, length(uniqs))
    }

    # For each query ratio...
    for (i in c(1:length(uniqs)))
    {
        # Query ratio
        ratio <- uniqs[[i]]
        
        # Reference ratio
        refRat <- refRats[[i]];
        
        t <- data[!is.na(data$ratio) & (data$ratio==ratio | data$ratio==refRat),]                

        # No FP or TP?
        if (length(unique(t$label)) == 1)
        {
            # No TP... Add a TP...
            if (unique(t$label) == 'FP')
            {
                t <- rbind(t, data.frame(label='TP', score=0, ratio=ratio))
            }
            
            # No FP... Add a FP...
            else
            {
                t <- rbind(t, data.frame(label='FP', score=0, ratio=ratio))
            }
        }
        
        t <- t[with(t, order(score)),]

        label <- ifelse(t$label == 'TP', 2, 1)
        preds <- prediction(t$score, label, label.ordering=c(1,2))
        perf  <- performance(preds, 'tpr', 'fpr')
        auc   <- performance(preds, 'auc')
        
        AUCs <- rbind(AUCs, data.frame(Ratio=ratio, AUC=round(unlist(auc@y.values), 4)))
        ROCs <- rbind(ROCs, data.frame(FPR=unlist(perf@x.values), TPR=unlist(perf@y.values), ratio=ratio))
    }

    ROCs$ratio = as.factor(ROCs$ratio)

    p <- ggplot(data=ROCs, aes_string(x='FPR', y='TPR'))             + 
            geom_abline(intercept=0, slope=1, linetype=2)            +
            geom_path(size=1, aes_string(colour='ratio'), alpha=0.5) +
            labs(colour=x$legTitle)                                  +
            theme_bw()

    if (!is.null(x$title))
    {
        p <- p + ggtitle(x$title)
    }
    
    if (!x$showLegend | length(uniqs) == 1)
    {
        p <- p + guides(colour=FALSE)
    }

    print(kable(AUCs))

    p <- .transformPlot(p)        
    print(p)

    return (list(AUC=AUCs))
}