\name{LadQuin}
\alias{LadQuin}
\title{LadQuin sequins}
\description{
Specify a design for ladder analysis.
}
\usage{
\method{LadQuin}{default}(mix)
}
\arguments{
  \item{seqs}{Unique sequin names}
  \item{input}{Input abundance (eg: attomol/ul)}
  \item{measured}{Measured abundance}
  \item{mix}{The name of the mixture file which the input concentration are to be read from. The file must be in the CSV format delimited by tabs. Each row of the file gives the input concentration of a sequin, which differs for each mixture type.

  The latest mixture files can be downloaded at www.sequin.xyz, it is important to use the correct mixture file. More details available at the website.
  
  Alternatively, \code{mix} can be 'A', 'B' or 'AB', for mixture A, B and AB respectively. Anaquin will download the mixture files, thus active Internet connection is required.}
}

\details{
  The \code{LadQuin} object combines the design information required for ladder analysis. While all arguments are optional, some information about the sequins must be given. Exactly what needs to be specified varies for by analysis, thus it is important to check the documentation carefully.

  Once created, the object can be used for ladder analysis. All ladder functions require this object.  
}
\value{
  An object of class \code{LadQuin} is returned. It contains the metadata required for data analysis. All ladder functions requires this object.
}

\author{Ted Wong \email{t.wong@garvan.org.au}}
\examples{
# Load a mixture from file
LadQuin(mix='mixture.csv')

# Load mixture A online
LadQuin(mix='A')

# Load mixture B online
LadQuin(mix='B')
}

\keyword{misc}
