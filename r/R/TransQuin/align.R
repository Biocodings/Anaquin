#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# ----------------------- Dilution Histogram -----------------------
#

plotDilutHist <- function()
{
    counts <- table(mtcars$gear)
    barplot(counts, main="",  xlab="Replicates", ylab = "Dilution %")
}

