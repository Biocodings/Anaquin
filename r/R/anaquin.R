#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

#
# Create a TransQuin data set for analyzing in Anaquin.
#

TransQuin <- function(mix=loadMixture(), ...)
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

    if (!is.null(x$baseMean))  { data$baseMean  <- x$baseMean  }  # TODO: Fix me
    if (!is.null(x$log2FoldChange)) { data$lfc <- x$log2FoldChange }  # TODO: Fix me
    if (!is.null(x$pvalue)) { data$pval <- x$pvalue }  # TODO: Fix me
    if (!is.null(x$expected.LFC)) { data$elfc <- x$expected.LFC }  # TODO: Fix me
    if (!is.null(x$lfcSE)) { data$lfcSE <- x$lfcSE }  # TODO: Fix me
    
    if (!is.null(x$X))  { data$X  <- x$X  }  # TODO: Fix me
    if (!is.null(x$A1)) { data$A1 <- x$A1 }  # TODO: Fix me
    if (!is.null(x$A2)) { data$A2 <- x$A2 }  # TODO: Fix me
    if (!is.null(x$A3)) { data$A3 <- x$A3 }  # TODO: Fix me
    if (!is.null(x$B1)) { data$B1 <- x$B1 }  # TODO: Fix me
    if (!is.null(x$B2)) { data$B2 <- x$B2 }  # TODO: Fix me
    if (!is.null(x$B3)) { data$B3 <- x$B3 }  # TODO: Fix me
    if (!is.null(x$prop)) { data$prop <- x$prop }  # TODO: Fix me

    r <- list('seqs'=data, mix=mix)
    class(r) <- 'TransQuin'

    return (r)
}

########################################################
#                                                      #
#                Accessor functions                    #
#                                                      #
########################################################

colors <- function(n)
{
    if (n == 5)
    {
        return (c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"))
    }
    else
    {
        return (c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd'))
    }
}

pval <- function(data)
{
    stopifnot(class(data) == 'TransQuin')
    
    if (is.null(data$seqs$pval))
    {
        # TODO: Implement me        
    }
    
    return (data$seqs$pval)
}

#
# Filtering functionality The following are supported:
#
#   'seqs'
#   'sequins'
#   'endos'
#
filter <- function(data, name)
{
    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin')
    
    stopifnot(name == 'seqs'    |
              name == 'sequins' |
              name == 'endos')
    
    if (name == 'seqs' | name == 'sequins')
    {
        if (!is.null(data$seqs$elf))
        {
            return (data$seqs[!is.na(data$seqs$elf),])
        }
        else
        {
            return (data$seqs) # TODO: Fix this
        }
    }
    else if (name == 'endos')
    {
        return (data$seqs[is.na(data$seqs$elf),])        
    }
}

sequins <- function(data)
{
    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin')

    return (row.names(filter(data, 'seqs')))
}

#
# Normalize the counts and calculate the base mean. Base mean is the average of the normalized count values, taken
# over all samples.
#

baseMean <- function(data)
{
    stopifnot(class(data) == 'TransQuin')
    
    if (is.null(data$seqs$baseMean))
    {
        # TODO: Implement me
    }
    
    return (data$seqs$baseMean)
}

names <- function(data)
{
    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin')

    return (row.names(data$seqs))
}

mLogFSE <- function(data, ids)
{
    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin')

    # Internal representation
    data <- data$seqs
    
    if (is.null(data$lfcSE))
    {
        # TODO: Implement me
    }
    
    return (data$lfcSE)
}

#
# Returns the measured logFold
#
mLogF <- function(data, ids)
{
    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin')
    
    if (is.null(data$seqs$lfc))
    {
        # TODO: Implement me
    }
    
    return (data$seqs$lfc)
}

#
# Returns the expected logFold. The following levels are supported:
#
#   TransQuin: 'exon'
#              'gene'
#              'isoform'
#

expectLF <- function(data, lvl, ids=NULL)
{
    stopifnot(!is.null(ids))

    stopifnot(lvl == 'gene'    ||
              lvl == 'isoform' ||
              lvl == 'exon')

    stopifnot(class(data) == 'TransQuin' ||
              class(data) == 'VarQuin'   ||
              class(data) == 'MetaQuin'  ||
              class(data) == 'Mixture')

    if (!is.null(data$seqs) & !is.null(data$seqs$elfc))
    {
        return (data$seqs$elf)
    }
    else
    {
        if (class(data) != 'Mixture')
        {
            data <- data$mix
        }
        
        if (is.null(data$elfc))
        {
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
        }
        
        return (r)
    }
}

#
# Returns the expected concentration for a sequin. The following levels are supported:
#
#   TransQuin: 'exon', 'isoform' and 'gene'
#

expectAbund <- function(data, ids, lvl, mix='A')
{
    stopifnot(lvl == 'exon'    |
              lvl == 'gene'    |
              lvl == 'isoform')

    stopifnot(class(data) == 'TransQuin' |
              class(data) == 'VarQuin'   |
              class(data) == 'MetaQuin'  |
              class(data) == 'Mixture')

    data <- data$mix
    
    switch(lvl, 'gene'    = { data <- data$genes[row.names(data$genes) %in% ids,]       },
                'isoform' = { data <- data$isoforms[row.names(data$isoforms) %in% ids,] },
                'exon'    = { data <- data$exons[row.names(data$exons) %in% ids,]       })

    if (is.null(data[mix]))
    {
        error(paste('Unknown mixture:', mix))
    }

    return (signif(data[mix][[1]], digits=2))
}

isoformsToGenes <- function(data, ids)
{
    stopifnot(class(data) == 'TransQuin')
  
    r <- data.frame(gID=rep(NA, length(ids))) 
    row.names(r) <- ids

    # Known isoforms
    knowns <- row.names(r)
    knowns <- knowns[row.names(r) %in% row.names(data$mix$isoforms)]

    # Map all the known isoforms    
    r[row.names(r) %in% knowns,] <- data$mix$isoforms[row.names(data$mix$isoforms) %in% knowns,]$GeneID

    return (r$gID)
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
        mix <- read.csv('/Users/tedwong/Sources/QA/data/trans/MTR004.v013.csv', row.names=1, sep='\t')
        #mix <- read.csv(url('https://s3.amazonaws.com/anaquin/mixtures/MTR004.v013.csv'), row.names=1, sep='\t')
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
