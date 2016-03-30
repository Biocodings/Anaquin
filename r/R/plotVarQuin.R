#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

plotAlleleAllele <- function(data)
{
    plotScatter(data, title='Expected allele frequency vs Measured allele frequency',
                      xname='Expected allele frequency',
                      yname='Measured allele frequency',
                  showStats='left',
                 showLegend=TRUE)
}

plotVAbundAbund <- function(data)
{
    plotScatter(data, title='Expected abundance vs Measured abundance',
                      xname='Expected abundance',
                      yname='Measured abundance',
                  showStats='left',
                 showLegend=FALSE)
}