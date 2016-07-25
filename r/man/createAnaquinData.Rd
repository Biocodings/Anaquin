\name{createAnaquinData}
\alias{createAnaquinData}
\title{Create Anaquin dataset}
\description{
Create an Anaquin dataset. The resulting R-object is required for every Anaquin R-function.
}
\usage{
\method{createAnaquinData}{default}(names, ...)
}
\arguments{
  \item{names}{Sequin names}
  \item{title}{Context specific inputs}
}
\details{
The function creates an R-object which every other Anaquin R-function requires. The object encapsulates the information required for analysis, thus the inputs are context specific. It is important to review the documentation for details. Users are encouraged to try our command-line Anaquin for auto-generating usable R-script.
}

\examples{
#
# Example 1: Create Anaquin dataset for PlotAssembly
#

data("Vignette_5.4.5")
data <- Vignette_5.4.5
CreateDataForAnaquin(names=row.names(data), input=data$input, sensitivity=data$Sn)
}

\author{Ted Wong \email{t.wong@garvan.org.au}}