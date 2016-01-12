#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

testExons <- function()
{
    m <- loadMixture()$exons

    #
    # The values are compared against Simon's computed fold-changes
    #
    
    checkTrue(m['R1_101:E001',]$fold == 0.125)
    checkTrue(m['R1_13:E008',]$fold  == 0.5)
    checkTrue(m['R1_42:E001',]$fold  == 0.125)
    checkTrue(m['R1_102:E008',]$fold == 0.00390625)
    checkTrue(m['R1_62:E007',]$fold  == 0.125)
    checkTrue(m['R1_63:E001',]$fold  == 4.0)
    checkTrue(m['R1_63:E002',]$fold  == 4.0)
}

testExons()