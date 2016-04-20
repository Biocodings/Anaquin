#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library('RUnit')
library('Anaquin')

testTransA <- function()
{
    x <- loadTransMix('A')
    
    checkEquals(class(x), 'TransMixture')
    checkEquals(nrow(x$genes), 78)
    checkEquals(nrow(x$isoforms), 164)
}

testTransAB <- function()
{
    x <- loadTransMix('AB')

    checkEquals(class(x), 'TransMixture')
    checkEquals(nrow(x$genes), 78)
    checkEquals(nrow(x$isoforms), 164)
}

testVarA <- function()
{
    x <- loadVarMix('A')
    
    checkEquals(class(x), 'VarMixture')
    checkEquals(nrow(x$seqs), 72)
}

testVarF <- function()
{
    x <- loadVarMix('F')
    
    checkEquals(class(x), 'VarMixture')
    checkEquals(nrow(x$seqs), 72)
}

#testTransAB()
#testVarA()
#testVarF()