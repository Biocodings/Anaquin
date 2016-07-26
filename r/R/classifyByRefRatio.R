#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

classifyByRefRatio <- function(inputs, refRatio)
{
    return (ifelse(inputs <= refRatio, 'FP', 'TP'))
}