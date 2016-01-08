#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# The point of inflection is defined as the level of concentration where accurate interpreation starts becoming questionable.
#

inflection <- function(x, y)
{
    require(ggplot2)
    
    #
    # For x=1,2,3,...,n, we would only fit b=3,4,...,n-2. Therefore the length of the frame is
    # n-4.
    #
    
    d <- data.frame(x=x, y=y)
    d <- d[order(x),]

    r <- data.frame(k=rep(NA,length(x)),
                    sums=rep(NA,length(x)),
                    lR2=rep(NA,length(x)),
                    rR2=rep(NA,length(x)))

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
        models <- plm(i)
        
        # Where this breakpoint occurs
        r$k[i] <<- models[[1]]
        
        m1 <- models[[2]]  # The left regression model
        m2 <- models[[3]]  # The right regression model

        r$lR2[i] <<- summary(m1)$r.squared
        r$rR2[i] <<- summary(m2)$r.squared

        # Calculate sum of residuals        
        r$sums[i] <<- sum((m1$residuals)^2) + sum((m2$residuals)^2)
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
    
    r <- list(m, r, b)
    r
}



