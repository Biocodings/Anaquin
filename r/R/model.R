#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Determine the soft detection limit for a given data set. The limit is defined as the level
# of concentration where accurate interpreation starts becoming questionable.
#

softLimit <- function(x, y)
{
    require(ggplot2)
    
    #
    # For x=1,2,3,...,n, we would only fit b=3,4,...,n-2. Therefore the length of the frame is
    # n-4.
    #
    
    d <- data.frame(x=x, y=y)
    d <- d[order(x),]
    r <- data.frame(k=rep(NA,length(x)), sums=rep(NA,length(x)))
    
    plm <- function(i)
    {
        #
        # Eg: (1,2,3,4) and i==2
        #
        #   -> (1,2) and (3,4). The index points to the last element in the first region as suggested
        #      in stats.stackexchange.com/questions/5700/finding-the-change-point-in-data-from-a-piecewise-linear-function
        #
        
        d1 <- head(d,i)
        d2 <- tail(d,-i)
        
        stopifnot(nrow(d1) >= 2)
        stopifnot(nrow(d2) >= 2)
        stopifnot(nrow(d1)+nrow(d2) == nrow(d))
        
        m1 <- lm(y~x, data=d1)
        m2 <- lm(y~x, data=d2)
        
        r <- list(d[i,]$x, m1, m2)
        r
    }
    
    lapply(2:(nrow(d)-3), function(i)
    {
        # Fit two piecewise linear models
        m <- plm(i)
        
        r$k[i]    <<- m[[1]]
        r$sums[i] <<- sum((m[[2]]$residuals)^2) + sum((m[[3]]$residuals)^2)
    })
    
    r <- r[!is.na(r$k),]
    
    #
    # Construct a plot of breaks vs total SSE
    #
    
    p <- ggplot(data = r, aes(x = k, y = sums))
    p <- p + xlab('Break point')
    p <- p + ylab('Total sum of squares')
    p <- p + geom_line()
    print(p)
    
    # The optimial breakpoint is where the minimum SSE is.
    b <- r[which.min(r$sums),]    
    
    # Fit the model again
    m <- plm(which.min(r$sums))
    
    r <- list(b, m, r)
    r
}
