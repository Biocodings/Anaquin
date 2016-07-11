#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.m2str <- function(m)
{
    eq <- substitute(italic(y) == a + b * italic(x)*','~~italic(r)^2~'='~r2, 
                     list(a  = format(coef(m)[1], digits = 2), 
                          b  = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

.lm2str <- function(data)
{
    return (.m2str(lm(y~x, data)))
}

#
# Limit-of-quantification (LOQ) is defined as the level of concentration where accurate interpreation starts becoming questionable.
#

showLOQ <- function(x, y, showDetails=FALSE)
{
    require(ggplot2)

    #
    # For x=1,2,3,...,n, we would only fit b=3,4,...,n-2. Therefore the length of the frame is n-4.
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
                    rInter=rep(NA,length(x)),
                    lr=rep(NA,length(x)),
                    rr=rep(NA,length(x)))

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
        
        return (list(breaks=d[i,]$x, 'lModel'=m1, 'rModel'=m2, 'lr'=cor(d1$x,d1$y), 'rr'=cor(d2$x,d2$y)))
    }

    lapply(2:(nrow(d)-3), function(i)
    {
        options(warn=-1)
        
        # Fit two piecewise linear models
        fit <- plm(i)
        
        options(warn=0)

        # Where this breakpoint occurs
        r$k[i] <<- fit$breaks

        m1 <- fit$lModel # The left regression model
        m2 <- fit$rModel # The right regression model

        sm1 <- summary(m1)
        sm2 <- summary(m2)
        
        r$lr[i]  <<- fit$lr
        r$rr[i]  <<- fit$rr
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

    #
    # How to define the LOQ breakpoint?
    #
    #   - Minimum where the SSE is
    #
    
    percentile <- ecdf(x)
    
    # Filter out the upper and lower quartiles
    a <- r[percentile(r$k) > 0.25 & percentile(r$k) <= 0.60,]
    
    # Where the minimum SSE is
    b1 <- r[which.min(r$sums),]    

    # Where the maximum above R2 is
    b2 <- a[which.max(a$rR2),]

    b1q <- percentile(b1$k) 
    b2q <- percentile(b2$k) 
    
    b <- b1

    if (b1q < 0.05)
    {
        # This is potentially a better breakpoint
        b <- b2        
    }

    # Fit the model again
    fit <- plm(as.numeric(row.names(b)))
    
    stopifnot(all.equal(summary(fit$lModel)$r.squared, b$lR2))  
    stopifnot(all.equal(summary(fit$rModel)$r.squared, b$rR2))
    stopifnot(all.equal(summary(fit$lModel)$coefficients[1,1], b$lInter))    
    stopifnot(all.equal(summary(fit$rModel)$coefficients[1,1], b$rInter))    
    stopifnot(all.equal(summary(fit$lModel)$coefficients[2,1], b$lSlope))    
    stopifnot(all.equal(summary(fit$rModel)$coefficients[2,1], b$rSlope))    
    
    data <- data.frame(x=x, y=y)

    #p <- ggplot(data=data, aes(x=x, y=y)) +
    #            geom_point() +
    #            geom_vline(xintercept=c(b$k), linetype="dotted") +
    #            theme_bw()
    #print(p)

    return (list('model'=fit, 'breaks'=b, 'details'=r))
}