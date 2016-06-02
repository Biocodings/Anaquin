\name{plotLOD}
\alias{plotLOD}
\title{Create LOD plot}
\description{
Create an Limit-of-Detection (LOD) plot.
}
\usage{
\method{plotLOD}{RNAQuin}(data, xlab, ylab, title='TransQuin Detection')

\method{plotLOD}{VarQuin}(data, xlab, ylab, title='VarQuin Detection')

\method{plotLOD}{MetaQuin}(data, xlab, ylab, title='MetaQuin Detection')

\method{plotLOD}{FusQuin}(data, xlab, ylab, title='FusQuin Detection')

\method{plotLOD}{LadQuin}(data, xlab, ylab, title='LadQuin Detection')
}
\arguments{
  \item{data}{Anaquin dataset. The data set needs to have input concentration and abundance.}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{title}{Title of the plot}
  \item{...}{For future expansion}  
}
\details{
Limit of detection plot indicates the confidence with which the features are detected relative to their expected abudnance (eg: input concentration).

For RNAQuin, the function emulates the LODR plot described in \emph{Munro 2014}:

------------------------------------------------------------------------------

\emph{
Identifying differentially expressed transcripts is the objective of differential expression experiments; however, how much information (signal) is needed to have confidence that a given fold change in expression of transcripts will be detected? With LODR estimates, empirical ERCC control ratio measurements can inform researchers of diagnostic power at a given fold change relative to transcript abundance for an experiment.

An LODR estimate for a particular fold change is the minimum signal, above which differentially expressed transcripts can be detected with a specified degree of confidence. LODR offers a statistically derived, objective alternative to other methods of parsing gene lists.

An LODR estimate is obtained for a specified ratio by modelling the relationship between differential expression test P-values and signal. An acceptable false discovery rate (FDR) must be chosen to estimate an LODR. For the selected FDR (q-value) a threshold P-value can be selected from the population of P-values from the experiment. An LODR estimate for each differential ratio is found based on the intersection of the model confidence interval upper bound with the P-value threshold.
}

------------------------------------------------------------------------------
}

\author{Ted Wong \email{t.wong@garvan.org.au}}
\examples{
#
# Example 1: RNAQuin
#
data(TransQuin)
plotLOD(data)
}