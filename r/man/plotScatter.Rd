\name{plotScatter}
\alias{plotScatter}
\title{Create scatter plot}
\description{
Create a scatter plot between expected abundance (eg: input concentration) against measured abundance (eg: FPKM).
}
\usage{
plotScatter(data, title, xlab, ylab, showSD, showLOQ, showAxis)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{createAnaquinData}. It needs to define information in \bold{Details}.}
  \item{title}{Label of the plot.}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{showSD}{Show standard deviation bars vertically? Default to FALSE.}
  \item{showLOQ}{Show limit-of-quantification? Default to TRUE.}
  \item{showAxis}{Show x-axis and y-axis? Default to TRUE.}  
}
\details{

\code{plotScatter} requires the following inputs from \code{createAnaquinData}.

\tabular{rlll}{
    \tab \code{names}       \tab Sequin names (eg: R1_1_1)\cr
    \tab \code{input/expected}       \tab Expected abundance. It must be either input or expected. input is preferred if the the x-axis is input concentration, while expected would be more comprehensible if the axis is not input concentration (eg: expected allele frequency).\cr
    \tab \code{measured} \tab Measured abundance (eg: FPKM, read counts)\cr
}

The scatter plot allows two variables be plotted; the x-axis is the expected abundance, and the y-axis is the measured abundance. Generally, the expected abundance is simply the input concentration, although other measures (eg: expected allele frequency) is also possible. The plot builds a linear regression between the two variables, and reports the standard statistics; R2, correlation and regression parameters on the graph.
}

\examples{
library("Anaquin")
data(seqCufflinks)

title <- 'Gene Expression'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'FPKM (log2)'

anaquin <- createAnaquinData(names=row.names(seqCufflinks), input=log2(seqCufflinks$InputConcent), measured=log2(seqCufflinks[,2:4]))
plotScatter(anaquin, title=title, xlab=xlab, ylab=ylab, showLOQ=TRUE)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}