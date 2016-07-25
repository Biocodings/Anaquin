#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test.isoformsToGenes <- function()
{
    seqs <- c('R1_1_1', 'R1_1_2', 'R1_1_3', 'R3_1_3')
    r <- .isoformsToGenes(seqs)

    checkEquals(length(r), 4)
    checkEquals(r[1], 'R1_1')
    checkEquals(r[2], 'R1_1')    
    checkEquals(r[3], 'R1_1')    
    checkEquals(r[4], 'R3_1')    
}

test.isoformsToGenes()