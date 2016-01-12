#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)

#
# Pool plot is a scatter plot where the variable in the x-axis is categorical. It is a common
# visualization tool for exporing the relationship sequin groups.
#

plotPool <- function(data,
                     shouldLog2 = TRUE,
                     xname = 'Expected log2 fold change of mixture A and B',
                     yname = 'Measured log2 fold change of mixture A and B')
{
    require(ggplot2)
    require(RColorBrewer)

    

}