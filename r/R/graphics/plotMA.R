#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute.
#

library(ggplot2)
library(gridExtra)
library(grid)

maData <- read.csv('/Users/tedwong/Sources/QA/r/maData.csv', row.names=1)
maDatAll <- read.csv('/Users/tedwong/Sources/QA/r/maDatAll.csv', row.names=1)
rm_dat <- read.csv('/Users/tedwong/Sources/QA/r/rm_dat.csv', row.names=1)
alphaPoint <- 0.8
ymalabel <- "Log2 Ratio of Normalized Counts"
xlabel = xlab("Log2 Average of Normalized Counts")
myXLimMA <- c(-14, 14)
myYLim <- c(-4, 4)

data <- read.csv('/Users/tedwong/Desktop/counts.txt', row.names=1)



#m <- read.csv('/Users/tedwong/Desktop/LODR_data_TED.csv', row.names=1)
#m <- m[!is.na(m$padj),]
#cutoffs <- plotLODR(row.names(m), m$baseMean, m$padj, m$expected.log2FoldChange)



ABCD <- function(data, mix=loadMixture())
{
    #
    # Bind the information for LODR with the data. We assume the status has already been classified.
    #

    data$status <- NA
    
    #
    # Eg:
    #      R1_43  2705.086868 above
    #      R1_52     7.499938 below
    #
    status <- cutoffs[[2]][c(1,4)]
    status <- status[!is.na(status$LODR),]
    
    for (i in 1:nrow(status))
    {
        t <- status[i,]
        
        if (nrow(data[data$Feature == row.names(t),]) > 0)
        {
            data[data$Feature == row.names(t),]$status <- t$LODR
        }
    }
    
    maStats <- function(x, c1, c2)
    {
        c(mean(log2(x[c1])-log2(x[c2])),  # M.Ave (used Y)
          sd(log2(x[c1])-log2(x[c2])),    # M.SD
          log2(mean(x)))                  # A     (used X)
    } 

    totCol <- 6 #ncol(data[-c(1:2)])
    
    #if(odd(totCol)) stop("Uneven number of replicates for the two sample types")

    maStatDat <- data.frame(t(apply(data[c(2:7)],
                                    1,
                                    maStats,
                                    c1 = c(1:(totCol/2)),
                                    c2 = c(((totCol/2)+1):totCol))))
    colnames(maStatDat) <- c("M.Ave","M.SD","A")
    
    data <- cbind(data, maStatDat, ratio=NA)
    data <- data[which(is.finite(data$M.Ave)),]    

    # Index for sequins
    si <- data$Feature %in% row.names(mix$genes)

    by(data[si,], 1:nrow(data[si,]), function(x)
    {
        logFold <- loadGene(x$Feature, mix)$logFold
        data[si,][data[si,]$Feature==x$Feature,]$ratio <<- abs(logFold)
    })
    
    data$ratio <- as.factor(as.character(data$ratio))

    # Now subset and continue with just sequins
    seqs <- data[si,]

    maPlot <- ggplot(seqs, aes(x = A, y = M.Ave)) +
               geom_point(data = subset(data, is.na(data$ratio)),
                          aes(x = A, y = M.Ave), colour = "grey80", alpha = 0.5) +
               
               geom_point(aes(colour = ratio), size = 5, alpha = alphaPoint) +
               geom_point(data = subset(data, status == 'below'), colour = "white", size = 2.5) +

               ylab(ymalabel) + xlabel +  coord_cartesian(xlim = c(-3,17), ylim = c(-10, 10)) +
               geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD, 
                          colour = ratio),size = 1,alpha = alphaPoint) +
               theme(legend.justification = c(1,0),legend.position=c(1,0)) +
               theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank()) +
               theme_bw() +
               scale_y_continuous(breaks = seq(-10, 10, 1))
    
#    tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
 #                        rowhead=list(fg_params = list(parse=TRUE)))

#    maPlot <- ggplot(maData, aes(x = A, y = M.Ave) ) + 
      #  geom_point(data = subset(maData, (LODR == "below")),
       #            colour = "white",size = 2.5) + 
    #    geom_hline(aes(yintercept = log2(Nominal), colour = Ratio), 
     #              alpha = 0.7) +
      #  geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
       #            size = 1, linetype = "longdash") + 

        #colScale + 
     #   annotation_custom(tableGrob(rm_dat,theme=tt), 
      #                    ymin = (myYLim[2]) - 0.25*myYLim[2], 
       #                   ymax = myYLim[2]) + 
    
    print(maPlot)
}

ABCD(data)







