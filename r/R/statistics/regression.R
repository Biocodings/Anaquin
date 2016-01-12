#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# The point of inflection is defined as the level of concentration where accurate interpreation starts becoming questionable.
#
#   - Piecewise linear
#

plotInflection <- function(x, y, method='piecewise', showDetails=FALSE)
{
    if (showDetails)
    {
        require(ggplot2)
    }

    #
    # For x=1,2,3,...,n, we would only fit b=3,4,...,n-2. Therefore the length of the frame is
    # n-4.
    #
    
    d <- data.frame(x=x, y=y)
    d <- d[order(x),]

    r <- data.frame(k=rep(NA,length(x)),
                    sums=rep(NA,length(x)),
                    lR2=rep(NA,length(x)),
                    lSlope=rep(NA,length(x)),                    
                    lInter=rep(NA,length(x)),
                    rR2=rep(NA,length(x)),
                    rSlope=rep(NA,length(x)),                    
                    rInter=rep(NA,length(x)))

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
        
        return (list(breaks=d[i,]$x, 'lModel'=m1, 'rModel'=m2))
    }

    lapply(2:(nrow(d)-3), function(i)
    {
        # Fit two piecewise linear models
        fit <- plm(i)
        
        # Where this breakpoint occurs
        r$k[i] <<- fit$breaks

        m1 <- fit$lModel  # The left regression model
        m2 <- fit$rModel  # The right regression model

        sm1 <- summary(m1)
        sm2 <- summary(m2)
        
        r$lR2[i] <<- sm1$r.squared
        r$rR2[i] <<- sm2$r.squared
        
        if (r$lR2[i] != 0)
        {
            r$lInter[i] <<- sm1$coefficients[1,1]
            r$lSlope[i] <<- sm1$coefficients[2,1]
        }
        
        if (r$rR2[i] != 0)
        {
            r$rInter[i] <<- sm2$coefficients[1,1]
            r$rSlope[i] <<- sm2$coefficients[2,1]
        }

        # Calculate sum of residuals        
        r$sums[i] <<- sum((m1$residuals)^2) + sum((m2$residuals)^2)
    })
    
    r <- r[!is.na(r$k),]
    
    #
    # Construct a plot of breaks vs total SSE
    #
    
    if (showDetails)
    {
        p <- ggplot(data = r, aes(x = k, y = sums))
        p <- p + xlab('Break point')
        p <- p + ylab('Total sum of squares')
        p <- p + geom_line()
        print(p)
    }

    # The optimial breakpoint is where the minimum SSE is.
    b <- r[which.min(r$sums),]    
    
    # Fit the model again
    fit <- plm(which.min(r$sums)+1)

    stopifnot(summary(fit$lModel)$r.squared == b$lR2)    
    stopifnot(summary(fit$rModel)$r.squared == b$rR2)
    stopifnot(summary(fit$lModel)$coefficients[1,1] == b$lInter)
    stopifnot(summary(fit$rModel)$coefficients[1,1] == b$rInter)    
    stopifnot(summary(fit$lModel)$coefficients[2,1] == b$lSlope)
    stopifnot(summary(fit$rModel)$coefficients[2,1] == b$rSlope)    
    
    return (list('model'=fit, 'breaks'=b, 'details'=r))
}
