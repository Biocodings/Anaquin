#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.plotROC <- function(data)
{
    require(ROCR)
    require(RColorBrewer)

    stopifnot(!is.null(data$pval) & !is.null(data$label) & !is.null(data$ratio))
    
    # Compute logarithm transformation to avoid overflowing
    data$lpval <- log2(data$pval)
    
    # Turn the probabilities into ranking classifer
    data$score <- 1-data$lpval
    
    ROCDat <- NULL
    AUCDat <- NULL
    
    # We'll render for each ratio
    ratios <- sort(data$ratio)
    
    for (ratio in unique(ratios))
    {
        #ratio <- 0.000244141
        
        if (!is.na(ratio))
        {
            t <- data[!is.na(data$ratio) & data$ratio == ratio,]
            
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
            
            label <- ifelse(t$label == 'TP', 2, 1)
            
            preds <- prediction(t$score, label, label.ordering=c(1,2))
            perf  <- performance(preds, 'tpr', 'fpr')
            auc   <- performance(preds, 'auc')
            
            # Now build the three vectors for plotting - TPR, FPR, and FoldChange
            AUC <- unlist(auc@y.values)
            
            print(paste(c('AUC for ', ratio, ': ', AUC), collapse = ''))
            
            AUCDatNew <- data.frame(ratio=ratio, AUC=round(AUC, digits=3))
            AUCDat <- rbind(AUCDat, AUCDatNew)
            
            FPR <- c(unlist(perf@x.values)) 
            TPR <- c(unlist(perf@y.values))
            
            ROCDatNew <- data.frame(FPR=FPR, TPR=TPR, ratio=ratio)
            ROCDat    <- rbind(ROCDat, ROCDatNew)
        }
    }
    
    ROCDat$ratio = as.factor(ROCDat$ratio)
    
    return (ROCDat)
}

.plotROC.Plot <- function(data, title=NULL)
{
    require(ggplot2)
    
    p <- ggplot(data=data, aes(x=FPR, y=TPR))                 + 
             geom_path(size=1, aes(colour=ratio), alpha=0.7)  + 
             geom_point(size=2, aes(colour=ratio), alpha=0.7) + 
             geom_abline(intercept=0, slope=1, linetype=2)    +
             labs(colour='Fold')                              +
             theme_bw()
    
    if (!is.null(title))
    {
        p <- p + ggtitle(title)
    }
    
    print(p)
}

plotROC.VarQuin <- function(data, title=NULL)
{
    .plotROC.Plot(.plotROC(data.frame(pval=data$seqs$pval, label=data$seqs$label, ratio=data$seqs$eAFreq)), title)
}

plotROC.TransQuin <- function(data, meth='validate')
{
    require(ROCR)
    require(ggplot2)
    
    stopifnot(meth == 'express' | meth == 'validate')
    stopifnot(class(data) == 'TransQuin' | class(data) == 'VarQuin')
    
    seqs  <- data$seqs
    lpval <- log2(seqs$pval)
    
    d <- data.frame(pval=seqs$pval, lpval=lpval, score=1-lpval, logFC=seqs$logFC)
    row.names(d) <- row.names(seqs)
    
    if (!is.null(seqs$cls))
    {
        d$cls <- seqs$cls    
    }
    
    ROCDat <- NULL
    AUCDat <- NULL
    
    logFCs <- d$logFC
    
    for (logFC in unique(logFCs))
    {
        if (logFC != 0)
        {
            t <- d[d$logFC == 0 | d$logFC == logFC,]
            
            if (meth == 'express')
            {
                stopifnot(!is.null(d$cls))
                
                t <- t[d$cls == 'TP' | d$cls == 'FP',]
                
                # 
                # From the reference manual:
                #
                # ... labels should be supplied as ordered factor(s), the lower level corresponding to the negative class, the upper level
                #     to the positive class ...
                #
                
                t$label <- ifelse(t$cls == 'TP', 2, 1)
                preds <- prediction(t$score, t$label, label.ordering=c(1,2))
            }
            
            #
            # Rank the p-values and group them like how's done for the ERCC dashboard. Refer to erccROC() in the ERCC dashboard
            # for reference implementation.
            #
            
            else
            {
                preds <- prediction(t$score, t$logFC, label.ordering=c(0, logFC))
            }
            
            perf  <- performance(preds, "tpr","fpr")
            auc   <- performance(preds, "auc")
            
            # Now build the three vectors for plotting - TPR, FPR, and FoldChange
            AUC <- unlist(auc@y.values)
            
            print(paste(c('AUC for ', logFC, ': ', AUC), collapse = ''))
            
            AUCDatNew <- data.frame(logFC=logFC, AUC=round(AUC, digits=3))
            AUCDat <- rbind(AUCDat, AUCDatNew)
            
            FPR <- c(unlist(perf@x.values)) 
            TPR <- c(unlist(perf@y.values))
            
            logFC <- c(rep(as.character(logFC), length(unlist(perf@y.values))))
            
            ROCDatNew <- data.frame(FPR=FPR, TPR=TPR, logFC=logFC)
            ROCDat    <- rbind(ROCDat, ROCDatNew)
        }
    }
    
    p <- ggplot(data=ROCDat, aes(x=FPR, y=TPR))              + 
            geom_path(size=2, aes(colour=logFC), alpha=0.7)  + 
            geom_point(size=5, aes(colour=logFC), alpha=0.7) + 
            geom_abline(intercept=0, slope=1, linetype=2)    +
            labs(colour='Log-Fold')                          +
            theme_bw()

    print(p)
}