\name{plotSensitivity}
\alias{plotSensitivity}
\title{Create sensitivity plot}
\description{
Create sensitivity plot where the x-axis is the logarithm of input concentration and the y-axis is the sensitivity.
}
\usage{
plotSensitivity(data, title, xlab, ylab, showLOA, threshold, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{createAnaquinData}. It needs to define information in \bold{Details}.}
  \item{title}{Label of the plot.}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{showLOA}{Show limit of assembly?}
  \item{threshold}{Threshold required for limit-of-assembly (LOA)}
  \item{...}{Reserved for internal testing}
}
\details{

\code{plotAssembly} requires the following inputs from \code{createAnaquinData}.

\tabular{rlll}{
    \tab \code{names}       \tab Sequin names\cr
    \tab \code{input}       \tab Input concentration (in attomole/ul)\cr
    \tab \code{sensitivity} \tab Measured sensitivity\cr
}

The sensitivity plot specifies the expected abundance (eg: input concentration) on the x-axis, and the measured sensitivity on the y-axis. Curve fitting is performed by non-linear least square fitting on sigmoid function. LOA is defined as the least abundant sequin reaching sensitivity of threshold.
}

\examples{
library("Anaquin")
data('seqCuffcompare')

title <- 'Assembly Plot'
xlab  <- 'Input Concentration (log2)'
ylab  <- 'Sensitivity'

names <- row.names(seqCuffcompare)
input <- log2(seqCuffcompare$InputConcent)
sensitivity <- seqCuffcompare$Sn

anaquin <- createAnaquinData(names=names, input=input, sensitivity=sensitivity)
plotSensitivity(anaquin, title=title, xlab=xlab, ylab=ylab, showLOA=TRUE)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}