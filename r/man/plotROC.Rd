\name{plotROC}
\alias{plotROC}
\title{Create ROC plot}
\description{
Create a receiver operating characteristic (ROC) plot. The curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at various threshold settings.
}
\usage{
\method{plotROC}{RNAQuin}(data, type=NULL, title='ROC Plot', color=NULL, guideTitle=NULL, showGuide=TRUE,  ...)

\method{plotROC}{VarQuin}(data, type=NULL, title='ROC Plot', color=NULL, guideTitle=NULL, showGuide=TRUE,  ...)

\method{plotROC}{MetaQuin}(data, type=NULL, title='ROC Plot', color=NULL, guideTitle=NULL, showGuide=TRUE,  ...)

\method{plotROC}{FusQuin}(data, type=NULL, title='ROC Plot', color=NULL, guideTitle=NULL, showGuide=TRUE,  ...)

\method{plotROC}{LadQuin}(data, type=NULL, title='ROC Plot', color=NULL, guideTitle=NULL, showGuide=TRUE,  ...)
}
\arguments{
  \item{data}{Anaquin data set. It needs to have input concentration and abundance.}}
  \item{type}{'SNP' or 'Indel' to filter the plot. Use \code{NULL} for all variants. }
  \item{title}{Title of the plot}
  \item{color}{Color of the curves, use \code{NULL} for default.}
  \item{guideTitle}{Title of the guide, use \code{NULL} for default.}
  \item{showGuide}{Show guides to the right?}
  \item{...}{For future expansion}
}
\author{Ted Wong \email{t.wong@garvan.org.au}}
\examples{
#
# Example 1: RNAQuin
#
data(RNAQuin)
plotROC(data)
}