#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

testVMixture <- function()
{
    mix <- loadMixture.VarQuin()

    checkEquals(0.5, alleleFreq(mix, 'D_1_1'))
}

testLFC <- function()
{
    mix <- loadMixture.TransQuin()
    r <- expectLF(mix, c('R2_7', 'R2_71', 'R2_72', 'R2_73', 'R2_76'), lvl='gene')

    checkEquals(0,  r[[1]])
    checkEquals(0,  r[[2]])
    checkEquals(-2, r[[3]])
    checkEquals(1,  r[[4]])
    checkEquals(3,  r[[5]])
}

testVMixture()
testLFC()