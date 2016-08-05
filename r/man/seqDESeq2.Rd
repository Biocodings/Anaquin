\docType{data}
\name{seqDESeq2}
\alias{seqDESeq2}
\title{Differential expression (RnaQuin)}
\format{A data frame:
\itemize{
  \item ExpLFC: Expected log-fold change 
  \item ObsLFC: Observed log-fold change
  \item SD: Standard deviation of the measurment
  \item Pval: P-value probability
  \item Qval: Q-value probability
  \item Mean: Average counts across the samples
  \item Label: Average counts across the samples  
}}
\usage{
data(seqDESeq2)
}
\description{
Differential gene expression estimated by DESeq2 from the RNA-Seq experiment of Hardwick (2016). The data correspond to the K562 and GM12878 samples.
}
\source{
S.A Hardwick. Spliced synthetic genes as internal controls in RNA sequencing experiments. \code{Nature Methods}, 2016.
}