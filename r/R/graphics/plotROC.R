#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotROC <- function(data, mix = loadMixture())
{
    require(ROCR)
    
    data <- data$seqs
    data <- data[data$class == 'TP' | data$class == 'FP',]
    
    # 
    # From the package's reference manual:
    #
    # ... labels should be supplied as ordered factor(s), the lower level corresponding to the negative class, the upper level
    #     to the positive class ...
    #

    data$scores <- 1 - data$qval
    data$label  <- ifelse(data$class == 'TP', 2, 1)
    
    pred <- prediction(data$scores, data$label, label.ordering=c(1,2))
    perf <- performance(pred, "tpr","fpr")
    
    plot(perf)
    
    return (list('pred' = pred, 'perf' = perf))
}
