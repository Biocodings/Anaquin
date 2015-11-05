#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Written by Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

#
# Determine the soft detection limit for a given data set. The limit is defined as the level
# of concentration where accurate interpreation starts becoming questionable.
#

softLimit <- function(x, y)
{
    #
    # Fit a piecewise linear model for each point. Note that this method always give the
    # optimial point, we should intrepret the result carefully. Generally speaking, it's
    # a proper soft limit if:
    #
    #   - Lowly expressed
    #   - The models are significant
    #   - The reduction in SSE is signficiant
    #   - Confirmed visually
    #
 
    #
    # For x=1,2,3,...,n, we would only fit b=2,3,...,n-1. Therefore the length of the frame is
    # n-2;
    #

    d <- data.frame(breaks=rep(0,length(x)-2),
                    models=rep(0,length(x)-2),
                    sums=rep(0,length(x)-2))
    
    lapply(2:(length(x)-2), function(i)
    {
        d$breaks <- i
    })    

    
    d
    
#    for (i in 1:(nrow(d)-2))
 #   {
  #      lm.shift(i+1)
  #  }
    
    
    
}