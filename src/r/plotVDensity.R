#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library(Anaquin)

bedgr <- '%3%/%4%'
annot <- '%3%/%5%'

data <- VarQuin(bedgr=bedgr, annot=annot)
plotCoverage(data)