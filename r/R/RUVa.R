#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#
#  Credits: RUV source package (http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)
#

#
# RNA normalization by RUVg (Remove Unwanted Variation). The code is slighlty different to the original
# R package. It includes positive and negative control, furthermore it also supports differential expression
# at the exon level.
#

.RUVaNorm <- function(x, cIdx, round=TRUE, k=1, epsilon=1, tolerance=1e-8, isLog=FALSE)
{
    # Log-linear GLM
    Y <- t(log(x+epsilon))
    
    # Scale to a zero-mean matrix
    Ycenter <- apply(Y, 2, function(x) scale(x, center = TRUE, scale=FALSE))
    
    m <- nrow(Y)
    n <- ncol(Y)
    
    # Perform a SVD decomposition for the control genes
    svdWa <- svd(Ycenter[, cIdx])
    
    drop <- 0
    first <- 1 + drop
    k <- min(k, max(which(svdWa$d > tolerance)))
    
    # Extract the first k principal components    
    W <- svdWa$u[, (first:k), drop = FALSE]
    
    # Estimate the parameter for the unwanted variation    
    alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
    
    # Estimate the biological variation
    correctedY <- Y - W %*% alpha
    
    if (!isLog)
    {
        if(round) {
            correctedY <- round(exp(correctedY) - epsilon)
            correctedY[correctedY<0] <- 0
        } else {
            correctedY <- exp(correctedY) - epsilon
        }
    }
    
    colnames(W) <- paste("W", seq(1, ncol(W)), sep="_")
    
    r <- list(W = W, normalizedCounts = t(correctedY))
    r
}
