library(RUnit)
library(Anaquin)

testPiecewiseLinear <- function()
{
    x <- c(-20:20)
    y <- abs(x)
    
    r <- plotInflection(x, y)

    checkEquals(r$breaks$k, -1)
}

testPiecewiseLinear()
