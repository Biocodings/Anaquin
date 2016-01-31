#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotROC <- function(data, meth='group')
{
    require(ROCR)
    require(ggplot2)
    
    stopifnot(meth == 'normal' | meth == 'group')
    stopifnot(class(data) == 'TransQuin' | class(data) == 'VarQuin')

    seqs <- data$seqs

    d <- data.frame(pvals=seqs$pval, scores=1-seqs$pval, logFC=seqs$logFC)
    
    ROCDat <- NULL
    AUCDat <- NULL
    
    if (meth == 'normal')
    {
        data <- data[data$cls == 'TP' | data$cls == 'FP',]
        
        # 
        # From the package's reference manual:
        #
        # ... labels should be supplied as ordered factor(s), the lower level corresponding to the negative class, the upper level
        #     to the positive class ...
        #
        
        data$scores <- 1 - pvals
        data$label  <- ifelse(data$cls == 'TP', 2, 1)
        
        pred <- prediction(data$scores, data$label, label.ordering=c(1,2))
        perf <- performance(pred, "tpr","fpr")
        
        plot(perf)
    }
    
    #
    # Rank the p-values and group them like how's done for the ERCC dashboard. Refer to erccROC() in the ERCC dashboard
    # for reference implementation.
    #
    
    else
    {
        for (logFC in unique(d$logFC))
        {
            if (logFC != 1)
            {
                t <- d[d$logFC == 1 | d$logFC == logFC,]
                
                preds <- prediction(t$scores, t$logFC, label.ordering=c(1, logFC))
                perf  <- performance(preds, "tpr","fpr")
                auc   <- performance(preds, "auc")
                
                # Now build the three vectors for plotting - TPR, FPR, and FoldChange
                AUC <- unlist(auc@y.values)
                
                AUCDatNew <- data.frame(logFC=logFC, AUC=round(AUC, digits=3))
                AUCDat <- rbind(AUCDat, AUCDatNew)
                
                FPR <- c(unlist(perf@x.values)) 
                TPR <- c(unlist(perf@y.values))
                
                logFC <- c(rep(as.character(logFC), length(unlist(perf@y.values))))
                
                ROCDatNew <- data.frame(FPR=FPR, TPR=TPR, logFC=logFC)
                ROCDat <- rbind(ROCDat, ROCDatNew)                
            }
        }
    }

    p <- ggplot(data=ROCDat, aes(x=FPR, y=TPR)) + 
         geom_path(size=2, aes(colour=logFC), alpha=0.7) + 
         geom_point(size=5, aes(colour=logFC), alpha=0.7) + 
        # colScale + 
        geom_abline(intercept=0, slope=1, linetype=2) +
         theme_bw() + 
       #  annotation_custom(grob=tableGrob(AUCAnnot, rows=NULL),
    #                      xmin=0.375, xmax=1.0, ymin=0, ymax=0.25) +
         theme(legend.position=c(0.75, 0.5))
    
    print(p)
    
    #return (list('pred' = pred, 'perf' = perf))
}