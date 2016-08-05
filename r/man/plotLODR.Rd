\name{plotLODR}
\alias{plotLODR}
\title{Create Limit-of-Detection Ratio (LODR) plot}
\description{
Create a Limit-of-Detection Ratio (LODR) plot between measured abundance (x-axis) and p-value probability (y-axis).
}
\usage{
plotLODR(data, title, xlab, ylab, FDR, legTitle, showConf)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{createAnaquinData}. It needs to define information in \bold{Details}.}
  \item{title}{Label of the plot. Default to NULL.}
  \item{xlab}{Label for the x-axis. Default to NULL.}
  \item{ylab}{Label for the y-axis. Default to NULL.}
  \item{FDR}{Chosen false-discovery-rate. Default to 0.1.}
  \item{legTitle}{Title for the legend. Default to ‘Ratio’.}
  \item{showConf}{Show confidence interval? Default to FALSE.}
}
\details{

\code{plotLODR} requires the following inputs from \code{createAnaquinData}.

\tabular{rlll}{
    \tab \code{names}    \tab Sequin names\cr
    \tab \code{measured} \tab Measured abundance (eg: average counts, DP field in a VCF file etc)\cr
    \tab \code{ratio}    \tab Expected ratio; eg: expected log-fold ratio or expected allele frequency etc\cr
    \tab \code{pval}     \tab P-value probability\cr
}

The LOD plot indicates the confidence in measuring synthetic features (eg: isoforms, SNPs, Indels, etc.) are detected relative to their abundance (eg: allele frequencies, read counts, etc.). The x-axis is the abundance, relative to the confidence (p-value) in the y-axis.
The curves are estimated by local regression estimations, and are colored by the sequin group, with 90% confidence interval indicated.
}

\examples{
library('Anaquin')
data('seqDESeq2')

xlab  <- 'Average Counts'
ylab  <- 'P-value'
title <- 'LODR Curves'

anaquin <- createAnaquinData(names=row.names(seqDESeq2), measured=seqDESeq2$Mean, ratio=seqDESeq2$ExpLFC, pval=seqDESeq2$Pval, qval=seqDESeq2$Qval)
plotLODR(anaquin, xlab=xlab, ylab=ylab, title=title, FDR=0.1)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}