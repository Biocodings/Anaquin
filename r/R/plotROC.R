#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotROC <- function(data, title=NULL, refRatio=NULL, shouldPseuoLog=TRUE, shouldAUC=FALSE)
{
    require(ROCR)
    require(grid)
    require(ggplot2)
    require(gridExtra)
    
    data <- data$seqs
    stopifnot(!is.null(data$pval) & !is.null(data$label))
    
    data <- data[, order(names(data))]
    
    # Compute logarithm transformation to avoid overflowing (also avoid pvalue of 0)
    if (shouldPseuoLog)
    {
        data$lpval <- log2(data$pval + 0.00001)
    }
    else
    {
        data$lpval <- log2(data$pval)
    }

    # Turn the probabilities into ranking classifer
    data$score <- 1-data$lpval

    if (is.null(data$ratio))
    {
        data$ratio <- abs(round(data$expected))
    }

    rocDat <- NULL
    aucDat <- NULL
    ratios <- sort(data$ratio)
    
    for (ratio in unique(ratios))
    {
        if (!is.na(ratio) && (is.null(refRatio) || ratio != refRatio))
        {
            if (is.null(refRatio))
            {
                t <- data[!is.na(data$ratio) & data$ratio == ratio,]
            }
            else
            {
                t <- data[!is.na(data$ratio) & (data$ratio == ratio | data$ratio == refRatio),]                
            }
            
            print(paste(c('Detectd for ', ratio, ': ', nrow(t)), collapse = ''))
            
            # No false-positive or true-positive?
            if (length(unique(t$label)) == 1)
            {
                # No TP... Add a TP...
                if (unique(t$label) == 'FP')
                {
                    x <- data.frame(label='TP', pval=0, ratio=ratio, lpval=0, score=0)
                }
            
                # No FP... Add a FP...
                else
                {
                    x <- data.frame(label='FP', pval=0, ratio=ratio, lpval=0, score=0)                    
                }
                
                t  <- rbind(t, x)
            }
            
            t <- t[with(t, order(score)),]

            #print(paste(c('Number of TP for ratio ', ratio, ':', nrow(t[t$label=='TP',])), collapse=' '))
            #print(paste(c('Number of FP for ratio ', ratio, ':', nrow(t[t$label=='FP',])), collapse=' '))

            label <- ifelse(t$label == 'TP', 2, 1)

            preds <- prediction(t$score, label, label.ordering=c(1,2))
            perf  <- performance(preds, 'tpr', 'fpr')
            auc   <- performance(preds, 'auc')

            AUC   <- round(unlist(auc@y.values), 4)
            print(paste(c('AUC for ', ratio, ': ', AUC), collapse=''))

            aucDatNew <- data.frame(Ratio=ratio, AUC=round(AUC, digits=3))
            aucDat <- rbind(aucDat, aucDatNew)
            
            FPR <- c(unlist(perf@x.values)) 
            TPR <- c(unlist(perf@y.values))
            
            rocDatNew <- data.frame(FPR=FPR, TPR=TPR, ratio=ratio)
            rocDat    <- rbind(rocDat, rocDatNew)
        }
    }
    
    rocDat$ratio = as.factor(rocDat$ratio)
    
    p <- ggplot(data=rocDat, aes(x=FPR, y=TPR))              + 
            geom_path(size=1, aes(colour=ratio), alpha=0.7)  + 
            geom_point(size=2, aes(colour=ratio), alpha=0.7) + 
            geom_abline(intercept=0, slope=1, linetype=2)    +
            labs(colour='Fold')                              +
            theme_bw()
    
    if (!is.null(title))
    {
        p <- p + ggtitle(title)
    }
    
    if (shouldAUC)
    {
        theme <- gridExtra::ttheme_default(core    = list(fg_params=list(cex = 0.7)),
                                           colhead = list(fg_params=list(cex = 1,0)),
                                           rowhead = list(fg_params=list(cex = 1.0)))
        g <- tableGrob(aucDat)
        p <- grid.arrange(p, g, ncol=1, heights=c(1.0,0.5))
    }

    print(p)    
}

plotROC.VarQuin <- function(data, title)
{
    .plotROC(data, title=title)
}

plotROC.TransQuin <- function(data, title)
{
    data$seqs <- TransDiff_(data)
    .plotROC(data, title=title, refRatio=0)
}

plotROC <- function(data, title=NULL)
{
    stopifnot (class(data) == 'TransQuin' || class(data) == 'VarQuin')
    
    if (class(data) == 'VarQuin')   { plotROC.VarQuin(data, title)   }
    if (class(data) == 'TransQuin') { plotROC.TransQuin(data, title) }
}