\name{VarQuin}
\alias{VarQuin}
\title{VarQuin sequins}
\description{
Specify a design for variant analysis.
}
\usage{
\method{VarQuin}{default}(mix)
}
\arguments{
  \item{mix}{the name of the mixture file which the input concentration are to be read from. The file must be in the CSV format delimited by tabs. Each row of the file gives the input concentration of a sequin, which differs for each mixture type.

  The latest mixture files can be downloaded at www.sequin.xyz, it is important to use the correct mixture file. More details available at the website.
  
  Alternatively, \code{mix} can be 'A', 'B' or 'AB', for mixture A, B and AB respectively. Anaquin will download the mixture files, thus active Internet connection is required.}

  \item{annot} {the name of the annotation file for variant analysis. The format of the file varies by analysis. Please check your analysis for further details.}
}

\details{
  The \code{VarQuin} object combines the design information required for variant analysis. While all arguments are optional, some information about the sequins must be given. Exactly what needs to be specified varies for by analysis, thus it is important to check the documentation carefully.

  Once created, the object can be used for variant analysis. All variant functions require this object.  
  
  The internal representation is unimportant and might change in future release as the command-line Anaquin stabilses its own API. 
}
\value{
  An object of class \code{FusQuin} is returned. It contains the metadata required for data analysis. All fusion functions requires this object.
}

\author{Ted Wong \email{t.wong@garvan.org.au}}
\note{The interface might change in future release as the command-line \code{Anaquin}
  stabilises its own API.}
\seealso
{
\code{\link{TransQuin}}

\code{\link{VarQuin}}
\code{\link{FusQuin}}
\code{\link{MetaQuin}}
}
\examples{
# Load a mixture from file
VarQuin(mix='/home/tedwong/mixture.csv')

# Load mixture A online
VarQuin(mix='A')

# Load mixture B online
VarQuin(mix='B')
}

\keyword{misc}
