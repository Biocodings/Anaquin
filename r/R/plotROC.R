#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotROC <- function(data, refRats, ...)
{
    require(ROCR)
    require(ggplot2)

    data <- data$seqs
    
    stopifnot(!is.null(refRats))
    stopifnot(!is.null(data$score))
    stopifnot(!is.null(data$label))
    stopifnot(!is.null(data$expected))    
    
    x <- list(...)
    
    if (is.null(x$showAUC))    { x$showAUC    <- FALSE   }
    if (is.null(x$legTitle))   { x$legTitle   <- 'Ratio' }
    if (is.null(x$showLegend)) { x$showLegend <- TRUE    }

    # This is the sequin groups
    data$ratio <- abs(round(data$expected))

    data <- data[data$label=='TP' | data$label=='FP',]
    data <- data[, order(names(data))]
    data <- data.frame(label=data$label, score=data$score, ratio=data$ratio)

    rocDat <- NULL
    aucDat <- NULL
    ratios <- sort(data$ratio)

    print('ROC Diagnostics: AUC')

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
        refRat <- ifelse(is.null(refRats), NULL, refRats[[i]]);
        
        if (is.null(refRat))
        {
            t <- data[!is.na(data$ratio) & data$ratio==ratio,]
        }
        else
        {
            t <- data[!is.na(data$ratio) & (data$ratio==ratio | data$ratio==refRat),]                
        }
        
        # No FP or TP?
        if (length(unique(t$label)) == 1)
        {
            # No TP... Add a TP...
            if (unique(t$label) == 'FP')
            {
                df <- data.frame(label='TP', score=0, ratio=ratio)
            }
            
            # No FP... Add a FP...
            else
            {
                df <- data.frame(label='FP', score=0, ratio=ratio)                    
            }
            
            t  <- rbind(t, df)
        }
        
        t <- t[with(t, order(score)),]

        label <- ifelse(t$label == 'TP', 2, 1)
        preds <- prediction(t$score, label, label.ordering=c(1,2))
        perf  <- performance(preds, 'tpr', 'fpr')
        auc   <- performance(preds, 'auc')
        
        AUC <- round(unlist(auc@y.values), 4)
        
        if (length(uniqs) == 1)
        {
            print(paste(c('AUC: ', AUC), collapse=''))
        }
        else
        {
            print(paste(c(ratio, ': ', AUC), collapse=''))
        }
        
        aucDatNew <- data.frame(Ratio=ratio, AUC=round(AUC, digits=3))
        aucDat <- rbind(aucDat, aucDatNew)
        
        FPR <- c(unlist(perf@x.values)) 
        TPR <- c(unlist(perf@y.values))
        
        rocDatNew <- data.frame(FPR=FPR, TPR=TPR, ratio=ratio)
        rocDat    <- rbind(rocDat, rocDatNew)
    }

    rocDat$ratio = as.factor(rocDat$ratio)
    
    p <- ggplot(data=rocDat, aes(x=FPR, y=TPR))           + 
            geom_abline(intercept=0, slope=1, linetype=2) +
            labs(colour=legTitle)                         +
            theme_bw()

    p <- p + geom_path(size=1, aes(colour=ratio), alpha=0.5)
    
    if (!is.null(title))
    {
        p <- p + ggtitle(title)
    }
    
    if (!x$showLegend | length(uniqs) == 1)
    {
        p <- p + guides(colour=FALSE)
    }

    p <- .transformPlot(p)        
    print(p)
}