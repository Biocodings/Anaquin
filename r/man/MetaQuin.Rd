\name{MetaQuin}
\alias{MetaQuin}
\title{MetaQuin sequins}
\description{
Specify a design for metagenomic analysis.
}
\usage{
\method{MetaQuin}{default}(mix)
}
\arguments{
  \item{mix}{the name of the mixture file which the input concentration are to be read from. The file must be in the CSV format delimited by tabs. Each row of the file gives the input concentration of a sequin, which differs for each mixture type.
  
  The latest mixture files can be downloaded at www.sequin.xyz, it is important to use the correct mixture file. More details at the website.
  
  Alternatively, \code{mix} can be 'A', 'B' or 'AB', for mixture A, B and AB respectively. Anaquin will download the mixture files, thus active Internet connection is required.}
}
\details{
  The \code{MetaQuin} object combines the design information required for metagenomic analysis. Theses objects are used by the statistical modelling and data analysis functions. While all arguments are optional, \code{mix} or \code{annot} must be given. If these variables are specified they must be a valid file either locally or remotely.
  
  Once created, the object can be used for metagenomic data analysis. All metagenomic functions will require this object.  
  
  The internal representation is unimportant and might change in future release as the command-line Anaquin stabilses its
  own API. 
}
\value{
  An object of class \code{MetaQuin} is returned. It contains the metadata required for data analysis. All metagenomic functions requires this object.
}

\author{Ted Wong \email{t.wong@garvan.org.au}}
\note{The interface might change in future release as the command-line \code{Anaquin}
  stabilises its own API.}
\examples{

# Load a mixture from file
MetaQuin(mix='/home/tedwong/mixture.csv')

# Load mixture A online
MetaQuin(mix='A')

# Load mixture B online
MetaQuin(mix='B')
}

\keyword{misc}
