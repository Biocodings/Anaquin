\name{FusQuin}
\alias{FusQuin}
\title{FusQuin sequins}
\description{
Specify a design for fusion analysis.
}
\usage{
\method{FusQuin}{default}(...)
}
\arguments{
  \item{seqs}{Unique sequin names (eg: F1_1_1)}
  \item{input}{Input abundance (eg: attomol/ul)}
  \item{measured}{Measured abundance (eg: number of spanning reads)}  
  \item{label}{Sequin classification}
  \item{mix}{the name of the mixture file which the input concentration are to be read from. The file must be in the CSV format delimited by tabs. Each row of the file gives the input concentration of a sequin, which differs for each mixture type.
  
  The latest mixture files can be downloaded at www.sequin.xyz, it is important to use the correct mixture file. More details available at the website.

\code{mix} can be 'A', for mixture A. Anaquin will download the mixture file, thus active Internet connection is required.}}

\details{
  The \code{FusQuin} object combines the design information required for fusion analysis. While all arguments are optional, some information about the sequins must be given. Exactly what needs to be specified varies by analysis, thus it is important to check the documentation carefully.

  Once created, the object can be used for fusion analysis. All fusion functions require this object.  
  The internal representation is unimportant and might change in future release. 
}
\value{
  An object of class \code{FusQuin} is returned. It contains the metadata required for data analysis. All fusion functions requires this object.
}

\author{Ted Wong \email{t.wong@garvan.org.au}}
\examples{
# Load a mixture from file
FusQuin(mix='mixture.csv')

# Load mixture A
FusQuin(mix='A')
}