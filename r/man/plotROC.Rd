\name{plotROC}
\alias{plotROC}
\title{Create ROC plot }
\description{
Create a receiver operating characteristic (ROC) plot. The curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at various threshold settings.
}
\usage{
\method{plotROC}{Anaquin}(data, refRats, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{CreateAnaquinData}. It needs to define information in \bold{Details}.}
  \item{title}{Label of the plot. Default to NULL.}
  \item{refRats}{Reference ratio group}
}
\details{

\code{plotROC} requires the following inputs from \code{CreateAnaquinData}.

\tabular{rlll}{
    \tab \code{names}       \tab Sequin names\cr
    
    \tab \code{label}       \tab Classified labels\cr
    \tab \code{score}       \tab How the ROC points should be ranked\cr    

    \tab \code{expected}       \tab Expected abundance; Input concentration (in attomole/ul).\cr
    \tab \code{measured}       \tab Measured abundance; how the sequins are ranked, for example, p-values and read depth etc.\cr
}

The ROC plot illustrates the performance of an experiment as the scoring threshold is varied. The curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at various threshold settings.

The plot provides diagnose performance to select possibly optimal models and to discard suboptimal ones. In particularly, the AUC statistics indicate the performance of the model relatively to a random experiment (AUC 0.5). 
}

\examples{

}

\author{Ted Wong \email{t.wong@garvan.org.au}}