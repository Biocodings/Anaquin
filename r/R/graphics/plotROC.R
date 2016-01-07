#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# ----------------------- ROC Plot -----------------------
#

plotROC <- function(r, labels, mix = loadMixture())
{
    require(ROCR)
    
    d <- r$data
    d <- d[!is.nan(d$measured),]
    d <- d[d$class == 'TP' | d$class == 'FP',]
    
    # 
    # From the package's reference manual:
    #
    # ... labels should be supplied as ordered factor(s), the lower level corresponding to the negative class, the upper level
    #     to the positive class ...
    #
    d$label <- ifelse(d$class == 'TP', 2, 1)
    
    d$scores <- 1 - d$pval
    
    require(ROCR)
    
    scores <- c(0.999985006663244,0.11848973121444,0.998317078130317,0.209432682723151,0.991844504940773,0.296175974698939,0.890247946913349,0.67291092334451,0.324408809166977,0.103365779196362,0.252794403922691,0.600822276080659,0.045599350395611,0.99997718730894,0.543686960607272,0.806223370592858,0.999698543215618,0.995835949956063,0.183529103779826,0.99986388387351,0.130236027425796,0.999552354071581,0.999917043997976,0.999911464673142,0.99985115274586,0.976245886879838,0.341646255148427,0.825269252328571,0.14267747036437,0.987840216112876,0.033142220590023)
    labels <- c(4,1,4,1,4,1,4,1,1,1,1,1,1,4,1,4,4,4,1,4,1,4,4,4,4,4,1,4,1,4,1)
    
    pred <- prediction(scores, labels, label.ordering=c(1,4))
    #pred <- prediction(d$scores, d$label, label.ordering=c(1,2))
    perf <- performance(pred, "tpr","fpr")
    
    #    plot(perf, colorize=TRUE, cex.lab=2, cex.main=2, lwd=10)    
    plot(perf)
}