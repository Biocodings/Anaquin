#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Create a TransQuin data set for analyzing in Anaquin.
#

transQuin <- function(mix=loadMixture(), ...)
{
    x <- list(...)
    
    # Sequin names must be present. We can use it to construct a data frame with the appropriate size.
    stopifnot(!is.null(x$seqs))
    
    data <- data.frame(row.names=x$seqs)
    
    if (!is.null(x$class))    { data$class    <- x$class    }
    if (!is.null(x$pval))     { data$pval     <- x$pval     }
    if (!is.null(x$qval))     { data$qval     <- x$qval     }
    if (!is.null(x$logFC))    { data$logFC    <- x$logFC    }    
    if (!is.null(x$ratio))    { data$ratio    <- x$ratio    }    
    if (!is.null(x$counts))   { data$counts   <- x$counts   }
    if (!is.null(x$expected)) { data$expected <- x$expected }
    if (!is.null(x$measured)) { data$measured <- x$measured }

    if (!is.null(x$X))  { data$X  <- x$X  }  # TODO: Fix me
    if (!is.null(x$A1)) { data$A1 <- x$A1 }  # TODO: Fix me
    if (!is.null(x$A2)) { data$A2 <- x$A2 }  # TODO: Fix me
    if (!is.null(x$A3)) { data$A3 <- x$A3 }  # TODO: Fix me
    if (!is.null(x$B1)) { data$B1 <- x$B1 }  # TODO: Fix me
    if (!is.null(x$B2)) { data$B2 <- x$B2 }  # TODO: Fix me
    if (!is.null(x$B3)) { data$B3 <- x$B3 }  # TODO: Fix me
    
    r <- list('seqs'=data, mix=mix)
    class(r) <- 'TransQuin'

    return (r)
}

########################################################
#                                                      #
#                Accessor functions                    #
#                                                      #
########################################################

#
# Returns the expected logFold. The following levels are supported:
#
#   TransQuin: 'exon'
#              'gene'
#              'isoform'
#

expectLF <- function(data, ids, lvl)
{
    stopifnot(lvl == 'gene'    |
              lvl == 'isoform' |
              lvl == 'exon')

    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin'  |
              class(data) == 'Mixture')

    if (class(data) != 'Mixture')
    {
        data <- data$mix
    }
    
    stopifnot(!is.null(data))

    #
    # Eg:
    #
    #      A        B      fold   logFold
    #   R2_59 0.4720688 0.4720688     1
    #

    switch(lvl, 'gene'    = { data <- data$genes    },
                'exon'    = { data <- data$exons    },
                'isoform' = { data <- data$isoforms })

    data <- data[row.names(data) %in% ids,]
    
    if (is.null(data$A) | is.null(data$B))
    {
        error(paste('Failed to find mixture A and B'))
    }

    r <- data.frame(logFC=round(log2(data$B / data$A)))
    row.names(r) <- row.names(data)
    return (r)
}

#
# Returns the expected concentration for a sequin. The following levels are supported:
#
#   TransQuin: 'exon', 'isoform' and 'gene'
#

expectAbund <- function(data, id, lvl, mix='A')
{
    stopifnot(lvl == 'gene' || lvl == 'isoform' || lvl == 'exon')
    stopifnot(class(data) == 'TransQuin' ||
              class(data) == 'VarQuin'   ||
              class(data) == 'MetaQuin'  ||
              class(data) == 'Mixture')

    data <- data$mix
    
    #
    # Eg:
    #
    #      A        B       fold   logFold
    #   R2_59 0.4720688 0.4720688     1
    #
    data <- data$genes[row.names(data$genes)==id,]

    stopifnot(nrow(data) <= 1)
        
    if (is.null(data[mix]))
    {
        error(paste('Unknown mixture:', mix))
    }

    return (signif(data[mix][[1]], digits=2))
}













.isoformsToGenes <- function(trans)
{
    trans <- as.character(trans)
    genes <- substr(as.character(trans), 1, nchar(trans)-2)
    genes
}

#
# Load the mixture into an R object that can be used in other Anaquin functions
#

loadMixture <- function(mix=NULL, exons=NULL)
{
    if (is.null(mix))
    {
        mix <- read.csv(url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv'), row.names=1, sep='\t')
        colnames(mix) <- c('Length', 'Mix.A', 'Mix.B')
    }

    if (is.null(exons))
    {
        #exons <- read.csv(url('https://s3.amazonaws.com/anaquin/mixtures/MTR006.v013.csv'), row.names=1)
        #exons$logFold <- as.numeric(as.character(exons$logFold))
    }
    
    mix$fold    <- mix$Mix.B / mix$Mix.A
    mix$logFold <- log(mix$fold)

    # Eg: R1_1 for R1_1_1 and R1_1_2    
    mix$GeneID <- .isoformsToGenes(row.names(mix))
    
    # Genes that are defined in the mixture
    geneIDs <- unique(mix$GeneID)
    
    g <- data.frame(row.names=geneIDs,
                    A=rep(0, length(geneIDs)),
                    B=rep(0, length(geneIDs)),
                    fold=rep(0, length(geneIDs)),
                    logFold=rep(0, length(geneIDs)),
                    type=rep(0, length(geneIDs)))
    
    #
    # Calculate the expected log-fold between mixture A and B
    #
    
    for (id in geneIDs)
    {
        seqs <- mix[mix$GeneID == id,]

        #
        # Calculate the expected abundance. We assume the following format:
        #
        #      ID     Length     Mix A      Mix B
        #    -------------------------------------
        #     R1_11    703     161.13281    5.0354
        #
        # 
        # We shouldn't assume anything for the column names.
        #

        g[id,]$A <- sum(seqs[,2])
        g[id,]$B <- sum(seqs[,3])
        
        #
        # Calculate the expected fold change (B/A)
        #
        
        g[id,]$fold    <- g[id,]$B / g[id,]$A
        g[id,]$logFold <- log2(g[id,]$fold)
        
        #
        # For each gene, we will classify it as:
        #
        #    'E' - equal concentration across mixtures
        #    'A' - concentration in the first mixture is higher
        #    'B' - concentration in the second mixture is higher
        #
        
        # We'll need that because we might have variable number of isoforms
        type = ''
        
        for (i in seqs$ID)
        {
            seq <- seqs[seqs$ID==i,]
            
            if (seq$Mix.A > seq$Mix.B)      { type <- paste(type, 'A', sep='') }
            else if (seq$Mix.B > seq$Mix.A) { type <- paste(type, 'B', sep='') }
            else                            { type <- paste(type, 'E', sep='') }
        }
        
        g[id,]$type <- type
    }
    
    i <- data.frame(mix)
    i <- i[with(i, order(row.names(i))),]

    r <- list('isoforms'=i, 'genes'=g[with(g, order(row.names(g))),], 'exons'=exons)
    class(r) <- c("Mixture")
    
    r$genes$logFold <- round(r$genes$logFold)
    r
}
