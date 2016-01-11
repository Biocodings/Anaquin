#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library(RUnit)
library(Anaquin)

test_1 <- function()
{
    d <- read.csv('tests/data/counts.txt', row.names=1)
    row.names(d) <- d$Feature
    d <- d[,-1]
    
    r <- plotMA(d)
    
    checkEquals(r$xname, 'Log2 Average of Normalized Counts')
    checkEquals(r$yname, 'Log2 Ratio of Normalized Counts')
}

test_1()