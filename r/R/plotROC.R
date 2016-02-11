#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotROC <- function(data, plotPerf=FALSE, refRatio=NULL, shouldPseuoLog=TRUE)
{
    require(ROCR)
    require(grid)
    require(gridExtra)
    
    stopifnot(!is.null(data$pval) & !is.null(data$label) & !is.null(data$ratio))
    
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
    
    ROCDat <- NULL
    AUCDat <- NULL
    
    # We'll render for each ratio
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
                # No TP...
                if (unique(t$label) == 'FP')
                {
                    x <- data.frame(pval=0, label='TP', ratio=ratio, lpval=0, score=0)
                }
            
                # No FP...
                else
                {
                    x <- data.frame(pval=0, label='FP', ratio=ratio, lpval=0, score=0)                    
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

            if (plotPerf)
            {
                tps <- unlist(perf@y.values)
                fps <- unlist(perf@x.values)

                # What's the FP when TP reaches 100%?
                cutoff_100 <- fps[which(tps == 1.0)[1]]
                
                # What's the FP when TP reaches 75%?
                cutoff_75  <- fps[which(tps >= 0.75)[1] - 1]                

                prec_75 = 0.75 / (0.75 + cutoff_75)
                
                print(paste(c('Precision for 75% for ratio', ratio, 'is :', prec_75), collapse=' '))
                
                plot(perf)
                mtext(paste(c('AUC:', AUC, 'for ratio:', ratio, '. TP==1.0, FP==', cutoff_100), collapse=' '))
            }

            AUCDatNew <- data.frame(Ratio=ratio, AUC=round(AUC, digits=3))
            AUCDat <- rbind(AUCDat, AUCDatNew)
            
            FPR <- c(unlist(perf@x.values)) 
            TPR <- c(unlist(perf@y.values))
            
            ROCDatNew <- data.frame(FPR=FPR, TPR=TPR, ratio=ratio)
            ROCDat    <- rbind(ROCDat, ROCDatNew)
        }
    }
    
    ROCDat$ratio = as.factor(ROCDat$ratio)
    
    return (list(roc=ROCDat, auc=AUCDat))
}

.plotROC.Plot <- function(data, title=NULL)
{
    require(ggplot2)
    
    aucData <- data$auc
    rocData <- data$roc
    
    p <- ggplot(data=rocData, aes(x=FPR, y=TPR))              + 
             geom_path(size=1, aes(colour=ratio), alpha=0.7)  + 
             geom_point(size=2, aes(colour=ratio), alpha=0.7) + 
             geom_abline(intercept=0, slope=1, linetype=2)    +
             labs(colour='Fold')                              +
             theme_bw()# + coord_fixed(ratio = 0.7) 
    
    if (!is.null(title))
    {
        p <- p + ggtitle(title)
    }

    g <- tableGrob(aucData)
    p <- grid.arrange(arrangeGrob(grobs = list(p, g)), ncol=1)
    
    print(p)
}

plotROC.VarQuin <- function(data, title=NULL, plotPerf=FALSE)
{
    .plotROC.Plot(.plotROC(data.frame(pval=data$seqs$pval, label=data$seqs$label, ratio=data$seqs$expected), plotPerf), title)
}

plotROC <- function(data, title=NULL, plotPerf=FALSE, refRatio=NULL, shouldPseuoLog=TRUE)
{
    stopifnot(class(data) == 'TransQuin' | class(data) == 'VarQuin')
    
    stopifnot(!is.null(data$seqs$pval))
    stopifnot(!is.null(data$seqs$label))
    stopifnot(!is.null(data$seqs$ratio))    
    stopifnot(!is.null(data$seqs$expected))

    data <- .plotROC(data.frame(pval=data$seqs$pval, label=data$seqs$label, ratio=data$seqs$expected), plotPerf=plotPerf,
                                                                                                       refRatio=refRatio,
                                                                                                 shouldPseuoLog=shouldPseuoLog)
    .plotROC.Plot(data, title)
}