\docType{data}
\name{seqCufflinks}
\alias{seqCufflinks}
\title{Gene expression (RnaQuin)}
\format{A data frame:
\itemize{
  \item InputConcent: Input concentration in attomol/ul
  \item Observed1: Measured FPKM for the first replicate
  \item Observed2: Measured FPKM for the second replicate
  \item Observed3: Measured FPKM for the third replicate  
}}
\usage{
data(seqCufflinks)
}
\description{
Gene expression estimated by Cufflinks from the RNA-Seq experiment of Hardwick (2016). The data correspond to the three replicate of the K562 cell lines.
}
\source{
S.A Hardwick. Spliced synthetic genes as internal controls in RNA sequencing experiments. \code{Nature Methods}, 2016.
}