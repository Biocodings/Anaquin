#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

testPiecewiseLinear <- function()
{
    x <- c(-20:20)
    y <- abs(x)
    
    r <- plotInflection(x, y)

    checkEquals(r$breaks$k, -1)
}

testPiecewiseLinear()
