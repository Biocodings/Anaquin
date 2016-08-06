\name{createAnaquinData}
\alias{createAnaquinData}
\title{Create Anaquin dataset}
\description{
Create an Anaquin dataset. The resulting R-object is required for every Anaquin R-function.
}
\usage{
createAnaquinData(...)
}
\arguments{
  \item{...}{Context specific inputs}
}
\details{
The function creates an R-object which every other Anaquin R-function requires. The object encapsulates the information required for analysis, thus the inputs are context specific. It is important to review the documentation for details. Users are encouraged to try our command-line Anaquin for auto-generating usable R-script.
}

\examples{
library("Anaquin")

#
# Example 1: Create Anaquin dataset for PlotSensitivity
#

data('seqCuffcompare')

title <- 'Assembly Plot'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'Sensitivity'

names <- row.names(seqCuffcompare)
input <- log2(seqCuffcompare$InputConcent)
sn <- seqCuffcompare$Sn

anaquin <- createAnaquinData(names=names, input=input, sensitivity=sn)
plotSensitivity(anaquin, title=title, xlab=xlab, ylab=ylab, showLOA=TRUE)

#
# Example 2: Create Anaquin dataset for PlotScatter
#

data(seqCufflinks)

title <- 'Gene Expression'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'FPKM (log2)'

names <- row.names(seqCufflinks)
input <- log2(seqCufflinks$InputConcent)
measured <- log2(seqCufflinks[,2:4])

anaquin <- createAnaquinData(names=names, input=input, measured=measured)
plotScatter(anaquin, title=title, xlab=xlab, ylab=ylab, showLOQ=TRUE)

#
# Example 3: Create Anaquin dataset for plotROC
#

data('seqDESeq2')

title <- 'ROC Plot'

names <- row.names(seqDESeq2)
ratio <- seqDESeq2$ExpLFC
measured <- seqDESeq2$ObsLFC
score <- 1-seqDESeq2$Pval
label <- seqDESeq2$Label

anaquin <- createAnaquinData(names=names, ratio=ratio, measured=measured, score=score, label=label)
plotROC(anaquin, title=title, refRats=0)

#
# Example 4: Create Anaquin dataset for plotLODR
#

data('seqDESeq2')

xlab  <- 'Average Counts'
ylab  <- 'P-value'
title <- 'LODR Curves'

names <- row.names(seqDESeq2)
measured <- seqDESeq2$Mean
ratio <- seqDESeq2$ExpLFC
pval <- seqDESeq2$Pval
qval <- seqDESeq2$Qval

anaquin <- createAnaquinData(names=names, measured=measured, ratio=ratio, pval=pval, qval=qval)
plotLODR(anaquin, xlab=xlab, ylab=ylab, title=title, FDR=0.1)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}