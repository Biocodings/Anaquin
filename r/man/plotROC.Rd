\name{plotROC}
\alias{plotROC}
\title{Create ROC plot }
\description{
Create a receiver operating characteristic (ROC) plot. The curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at various threshold settings.
}
\usage{
plotROC(data, refRats, title, legTitle, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{createAnaquinData}. It needs to define information in \bold{Details}.}
  \item{title}{Label of the plot. Default to NULL.}
  \item{legTitle}{Title of the legend. Default to NULL.}
  \item{refRats}{Reference ratio groups. Default to 'Ratio'.}
  \item{...}{Reserved for internal testing}
}
\details{

\code{plotROC} requires the following inputs from \code{createAnaquinData}.

\tabular{rlll}{
    \tab \code{names}       \tab Sequin names\cr
    \tab \code{label}       \tab Classified labels\cr
    \tab \code{score}       \tab How the ROC points should be ranked\cr    
    \tab \code{ratio}       \tab Expected ratio; eg: expected log-fold ratio\cr
}

The ROC plot illustrates the performance of an experiment as the scoring threshold is varied. The curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at various threshold settings.

The plot provides diagnose performance to select possibly optimal models and to discard suboptimal ones. In particularly, the AUC statistics indicate the performance of the model relatively to a random experiment (AUC 0.5). 
}

\examples{
library('Anaquin')
data('seqDESeq2')

title <- 'ROC Plot'

# Sequin names
names <- row.names(seqDESeq2)

# Expected log-fold
ratio <- seqDESeq2$ExpLFC

# Measured log-fold
measured <- seqDESeq2$ObsLFC

# How the ROC curves are ranked
score <- 1-seqDESeq2$Pval

# Classified labels
label <- seqDESeq2$Label

anaquin <- createAnaquinData(names=names, ratio=ratio, score=score, label=label)
plotROC(anaquin, title=title, refRats=0)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}