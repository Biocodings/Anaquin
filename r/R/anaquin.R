#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.createData <- function(x, keys)
{
    # Do we have sequin names?
    if (is.null(x$names))
    {
        return (NULL)
    }
    
    data <- data.frame(row.names=x$names)
    
    for (key in keys)
    {
        if (!is.null(x[[key]]))
        {
            data[[key]] <- x[[key]]
        }
    }

    return (data)    
}

CreateDataForAnaquin <- function(...)
{
    x <- list(...)

    # Metrics Anaquin supports
    keys <- c('sd', 'pval', 'qval', 'ratio', 'mean', 'expected', 'input', 'measured', 'sensitivity', 'label', 'score')

    r <- list('seqs'=.createData(x, keys))
    
    if (is.null(r$seqs))
    {
        stop('Please specify sequin names, for example, CreateDataForAnaquin(names=...)')
    }
    
    class(r) <- 'Anaquin'
    return (r)
}

#
# Convert RnaQuin sequin isoform to sequin gene. The inputs assumed be valid.
#
#    Eg: R1_1_1 to R1_1
#
.isoformsToGenes <- function(ids)
{
    ids <- strsplit(ids, '_')
    
    f <- function(x)
    {
        paste(x[1:length(x)-1], collapse='_')
    }
    
    unlist(lapply(ids, f))
}