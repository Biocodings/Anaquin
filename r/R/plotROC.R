#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotROC <- function(data,
                    legTitle='LFC',
                    title=NULL,
                    refRats=NULL,
                    showGuide=TRUE,
                    ...)
{
    require(ROCR)
    require(grid)
    require(ggplot2)

    data <- data$seqs

    stopifnot(!is.null(data$score))
    stopifnot(!is.null(data$label))
    
    if (is.null(data$expected))
    {
        data$expected <- 1
    }
    
    showAUC <- FALSE
    
    if (is.null(data$ratio))
    {
        data$ratio <- abs(round(data$expected))
    }

    data <- data[data$label=='TP' | data$label=='FP',]
    data <- data[, order(names(data))]
    data <- data.frame(label=data$label, score=data$score, ratio=data$ratio)

    rocDat <- NULL
    aucDat <- NULL
    ratios <- sort(data$ratio)

    print('ROC Diagnostics: AUC')

    uniqs <- unique(ratios)
    uniqs <- uniqs[!(uniqs %in% refRats)]

    #
    # The reference ratios need to have the same length as the query ratios. We can do some sanity
    # checks to match the size.
    #
    
    if (length(refRats) == 1 && length(refRats) != length(uniqs))
    {
        refRats <- rep(refRats, length(uniqs))
    }
    
    # For each ratio...
    for (i in c(1:length(uniqs)))
    {
        # Query ratio
        ratio <- uniqs[[i]]
        
        # Reference ratio
        refRat <- ifelse(is.null(refRats), NULL, refRats[[i]]);
        
        if (is.null(refRat))
        {
            t <- data[!is.na(data$ratio) & data$ratio== ratio,]
        }
        else
        {
            t <- data[!is.na(data$ratio) & (data$ratio == ratio | data$ratio == refRat),]                
        }
        
        # No FP or TP?
        if (length(unique(t$label)) == 1)
        {
            # No TP... Add a TP...
            if (unique(t$label) == 'FP')
            {
                x <- data.frame(label='TP', score=0, ratio=ratio)
            }
            
            # No FP... Add a FP...
            else
            {
                x <- data.frame(label='FP', score=0, ratio=ratio)                    
            }
            
            t  <- rbind(t, x)
        }
        
        t <- t[with(t, order(score)),]
        #t[t$label=='FP',]$score <- -1000
        
        #print(paste(c('Number of TP for ratio ', ratio, ':', nrow(t[t$label=='TP',])), collapse=' '))
        #print(paste(c('Number of FP for ratio ', ratio, ':', nrow(t[t$label=='FP',])), collapse=' '))
        
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
    
    if (showAUC)
    {
        theme <- gridExtra::ttheme_default(core    = list(fg_params=list(cex = 0.7)),
                                           colhead = list(fg_params=list(cex = 1,0)),
                                           rowhead = list(fg_params=list(cex = 1.0)))
        g <- tableGrob(aucDat)
        p <- grid.arrange(p, g, ncol=1, heights=c(1.0,0.5))
    }

    if (!showGuide | length(uniqs) == 1)
    {
        p <- p + guides(colour=FALSE)
    }

    p <- .transformPlot(p)        
    print(p)
}

#plotROC.FusQuin <- function(data, title, color, type)
#{
#    if (is.null(title)) { title <- 'FusQuin Detection' }
#    
#    data$seqs$pval <- (max(data$seqs$measured) + 1) - data$seqs$measured
#    .plotROC(data, title=title, color=color, refRats=0, showGuide=FALSE)
#}

#plotROC.TransQuin <- function(data, title, color, type)
#{
#    data$seqs <- TransDiff_(data)
#    .plotROC(data, title=title, color=color, refRats=0)
#}