#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Written by Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#
#  Credits: RUV source package (http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)
#

#
# Apply RNA normalization to a count matrix. The following modes are supported:
#
#    'RUV' -> Remove Unwanted Variation (http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)
#

.isWholeNumber <- function(x, tol = .Machine$double.eps^0.5)
{
    abs(x - round(x)) < tol
}

TransNorm <- function(r, m, mode, k=1)
{
    if (!.isWholeNumber(x))
    {
        stop(paste0("The count matrix should contain only positive numbers."))
    }
    
    Y <- t(x)

    Ycenter <- apply(Y, 2, function(x) scale(x, center = TRUE, scale=FALSE))

    m <- nrow(Y)
    n <- ncol(Y)

    svdWa <- svd(Ycenter[, r$genes])
    
    drop <- 0
    
    first <- 1 + drop
    k <- min(k, max(which(svdWa$d > tolerance)))
    W <- svdWa$u[, (first:k), drop = FALSE]
    alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
    correctedY <- Y - W %*% alpha
    
    if (!isLog & all(.isWholeNumber(x)))
    {
        if(round) {
            correctedY <- round(exp(correctedY) - epsilon)
            correctedY[correctedY<0] <- 0
        } else {
            correctedY <- exp(correctedY) - epsilon
        }
    }

    colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
    return(list(W = W, normalizedCounts = t(correctedY)))
}