#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotTROC <- function(data, title=NULL, plotPerf=FALSE, refRatio=NULL, shouldPseuoLog=TRUE)
{
    # Classify the sequins
    data$seqs <- TransDiff(data)

    plotROC(data, title=title, plotPerf=plotPerf, refRatio=0, shouldPseuoLog=shouldPseuoLog)
}

plotTLogFold <- function(data)
{
    plotScatter(data, title='Expected log-fold vs Measured log-fold',
                xname='Expected log-fold',
                yname='Measured log-fold',
                shouldLog2=FALSE,
                showStats='left',
                showLegend=FALSE)
}
