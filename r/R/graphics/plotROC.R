#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# For a description of the ROC curve, try: http://ebp.uga.edu/courses/Chapter%204%20-%20Diagnosis%20I/8%20-%20ROC%20curves.html
#

plotROC <- function(data, mix = loadMixture())
{
    require(ROCR)

    stopifnot(class(data) == 'TransQuin')

    data <- data$seqs
    
    # The input is expected to be classifed (eg: TransClassify)
    stopifnot(is.null(data$seqs))

    data <- data[data$class == 'TP' | data$class == 'FP',]
    
    # 
    # From the package's reference manual:
    #
    # ... labels should be supplied as ordered factor(s), the lower level corresponding to the negative class, the upper level
    #     to the positive class ...
    #

    data$scores <- 1 - data$pval
    data$label  <- ifelse(data$class == 'TP', 2, 1)

    pred <- prediction(data$scores, data$label, label.ordering=c(1,2))
    perf <- performance(pred, "tpr","fpr")
    
    plot(perf)
    
    return (list('pred' = pred, 'perf' = perf))
}
