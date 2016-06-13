#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Classification for differential analysis. Classify each differential test by:
#
#    TP: fold-change with more than logFC and expressed
#    FP: fold-change with at most logFC and expressed
#    FN: fold-change with more than logFC and not expressed
#    TN: fold-change with at most logFC and not expressed
#
TransDiff_ <- function(data, qCutoff=0.1, logFC=0)
{
    print(paste('Probability threshold:', qCutoff))

    data <- data$seqs
    data$label <- NA

    stopifnot(!is.null(data$expected))
    
    # Expected log-fold
    data$elfc <- data$expected #expectLF(data, lvl=lvl, ids=row.names(seqs))$logFC
    
    for (id in row.names(data))
    {
        qval = data[id,]$qval

        if  (!is.na(qval))
        {
            elfc <- data[id,]$elfc
            
            if (length(elfc) > 0)
            {
                #
                # Say if the known log-fold change is -3, should this be differentially expressed? That depends on the
                # threshold. Usuaully, we'd assume anything more than log-fold change of 0 would be expressed.
                #
                
                # Differential expressed?
                #if (qval <= qCutoff)
                {
                    data[id,]$label <- ifelse(abs(elfc) <= abs(logFC), 'FP', 'TP')
                }

                # Non-dfifferential expressed?
                #else
                #{
                #    data[id,]$label <- ifelse(abs(elfc) <= abs(logFC), 'TN', 'FN')
                #}
            }
        }
    }
    
    print(sprintf("Detected %d false positives", nrow(data[data$label=='FP',])))
    print(sprintf("Detected %d true positives",  nrow(data[data$label=='TP',])))
    print(sprintf("Detected %d true negatives",  nrow(data[data$label=='TN',])))
    print(sprintf("Detected %d false negatives", nrow(data[data$label=='FN',])))

    return (data)
}

.TransDiff <- function(data, r, detected, p, logFC)
{
    haveMean <- !is.null(r$mean)
    
    # Reference genes    
    refs <- as.character(row.names(data$mix$genes))
    
    # Genes that have been detected in the experiment
    detected <- rownames(r) %in% refs
    
    # Filter out to only the rows with sequins
    r <- r[detected,]
    
    print(sprintf("Detected %d reference sequins", length(refs)))
    print(sprintf("Detected %d experimental genes", length(rownames(r))))
    print(sprintf("%d sequins failed to detect", length(refs) - length(rownames(r))))
    
    #
    # Create a data-frame for all sequins, whether it's detected. The fold-changes are
    # on the logarithmic scale.
    #
    # NaN refers to undetected sequins while NA refers to detected but untested sequins
    #
    
    n <- rep(NaN, length(refs))
    x <- data.frame(mean=n,
                    expected=n,
                    measured=n,
                    pval=n,
                    qval=n,                    
                    expressed=n,
                    class=n)
    rownames(x) <- refs
    
    mix <- data$mix
    
    for (id in refs)
    {
        x[id,]$expected <- mix$genes[row.names(mix$genes)==id,]$logFold
    }
    
    for (id in rownames(r))
    {
        if (!is.null(r[id,]$mean))
        {
            x[id,]$mean <- r[id,]$mean
        }

        x[id,]$pval      <- r[id,]$pval
        x[id,]$qval      <- r[id,]$qval
        x[id,]$measured  <- r[id,]$lfc
        x[id,]$expressed <- ifelse(r[id,]$qval <= p, 'T', 'F')
    }

    # Sort by adjusted p-values so that the expressed genes are at the front
    x <- x[with(x, order(qval)),]

    # Fit a simple-linear regression model
    m <- lm(measured ~ expected, x[is.finite(x$measured),])
    
    r     <- cor(as.numeric(x$expected), as.numeric(x$measured))
    r2    <- summary(m)$r.squared
    slope <- coef(m)["refs"]
    
    #d <- .classify(d, logFC)

    x  <- x[!is.na(x$pval),]
    x  <- x[!is.nan(x$pval),]
    td <- TransQuin(seqs=row.names(x), expected=x$expect, measured=x$measured, pval=x$pval, qval=x$qval)

    if (haveMean)
    {
        td$seqs$mean <- x$mean
    }

    plotExpress(td,
                showLOQ=FALSE,
                title='TransQuin Differential',
                xlab='Expected log-fold (log2)',
                ylab='Measured log-fold (log2)')

    plotROC(td, title='TransQuin Differential')
    
    if (haveMean)
    {
        plotLODR(td, chosenFDR=0.1)
    }
}

TransDiff.DESeq2 <- function(data, r, p=0.1, logFC=0)
{
    require(DESeq2)

    stopifnot(class(data) == 'TransQuin')
    
    r <- data.frame(mean=r$baseMean,
                    lfc=r$log2FoldChange,
                    pval=r$pvalue,
                    qval=r$padj,
                    row.names=rownames(r))

    .TransDiff(data, r, rownames(r), p, logFC)
}

TransDiff.edgeR <- function(data, r, p=0.1, logFC=0)
{
    require(edgeR)
    
    stopifnot(class(data) == 'TransQuin')
    
    r <- data.frame(lfc=r$table$logFC,
                    pval=r$table$PValue,
                    qval=r$table$PValue,
                    row.names=rownames(r))

    .TransDiff(data, r, rownames(r), p, logFC)    
}

TransDiff <- function(data, r)
{
    stopifnot(class(r) == 'DESeqResults' || class(r) == 'DGEExact')
    
    if (class(r) == 'DESeqResults') { TransDiff.DESeq2(data, r) }
    if (class(r) == 'DGEExact')     { TransDiff.edgeR(data, r)  }
}
