\name{plotSensitivity}
\alias{plotSensitivity}
\title{Create sensitivity plot}
\description{
Create sensitivity plot where the x-axis is the logarithm of input concentration and the y-axis is the sensitivity.
}
\usage{
\method{plotSensitivity}{Anaquin}(data, ...)
}
\arguments{
  \item{data}{Anaquin dataset created by \code{CreateAnaquinData}. It needs to define information in \bold{Details}.}
  \item{title}{Label of the plot.}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{showLOA}{Show limit of assembly?}
  \item{threshold}{Threshold required for limit-of-assembly (LOA)}
}
\details{

\code{plotAssembly} requires the following inputs from \code{CreateAnaquinData}.

\tabular{rlll}{
    \tab \code{names}       \tab Sequin names\cr
    \tab \code{input}       \tab Input concentration (in attomole/ul)\cr
    \tab \code{sensitivity} \tab Measured sensitivity\cr
}

The sensitivity plot specifies the expected abundance (eg: input concentration) on the x-axis, and the measured sensitivity on the y-axis. Curve fitting is performed by non-linear least square fitting on sigmoid function. LOA is defined as the least abundant sequin reaching sensitivity of threshold.
}

\examples{
library("Anaquin")

data("Vignette_5.4.5")
data <- Vignette_5.4.5
data <- createAnaquinData(names=row.names(data), input=log2(data$InputConcent), sensitivity=data$Sn)

plotSensitivity(data)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}