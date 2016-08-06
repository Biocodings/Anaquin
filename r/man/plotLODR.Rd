\name{plotLODR}
\alias{plotLODR}
\title{Create Limit-of-Detection Ratio (LODR) plot}
\description{
Create a Limit-of-Detection Ratio (LODR) plot between measured abundance (x-axis) and p-value probability (y-axis).
}
\usage{
plotLODR(data, FDR, title, xlab, ylab, legTitle, showConf, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{createAnaquinData}. It needs to define information in \bold{Details}.}
  \item{FDR}{Chosen false-discovery-rate. Default to \code{0.1}.}
  \item{title}{Label of the plot. Default to \code{NULL}.}
  \item{xlab}{Label for the x-axis. Default to \code{NULL}.}
  \item{ylab}{Label for the y-axis. Default to \code{NULL}.}
  \item{legTitle}{Title for the legend. Default to \code{'Ratio'}.}
  \item{showConf}{Show confidence interval? Default to \code{FALSE}.}
  \item{...}{Reserved for internal testing}
}
\details{

\code{plotLODR} requires the following inputs from \code{createAnaquinData}.

\tabular{rlll}{
    \tab \code{names}    \tab Sequin names\cr
    \tab \code{measured} \tab Measured abundance (eg: average counts, DP field in a VCF file etc)\cr
    \tab \code{ratio}    \tab Expected ratio; eg: expected log-fold ratio\cr
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

# Sequin names
names <- row.names(seqDESeq2)

# Expected log-fold
ratio <- seqDESeq2$ExpLFC

# Measured average abundance
measured <- seqDESeq2$Mean

pval <- seqDESeq2$Pval
qval <- seqDESeq2$Qval

anaquin <- createAnaquinData(names=names, measured=measured, ratio=ratio, pval=pval, qval=qval)
plotLODR(anaquin, xlab=xlab, ylab=ylab, title=title, FDR=0.1)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}