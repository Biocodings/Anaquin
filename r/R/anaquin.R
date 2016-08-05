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

createAnaquinData <- function(...)
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
isoformsToGenes <- function(names)
{
    names <- strsplit(names, '_')
    
    f <- function(x)
    {
        paste(x[1:length(x)-1], collapse='_')
    }
    
    unlist(lapply(names, f))
}


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