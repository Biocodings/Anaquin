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


tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     rowhead=list(fg_params = list(parse=TRUE)))


maPlot <- ggplot(maData, aes(x = A, y = M.Ave) ) + 
    geom_point(data = subset(maDatAll, maDatAll$Ratio == "Endo"),
               aes(x = A, y = M.Ave),
               colour = "grey80", alpha = 0.5) +
    geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD, 
                      colour = Ratio),size = 1,alpha = alphaPoint) + 
    geom_point(aes(colour = Ratio),size = 5, alpha = alphaPoint) +
    geom_point(data = subset(maData, (LODR == "below")),
               colour = "white",size = 2.5) + 
    geom_hline(aes(yintercept = log2(Nominal), colour = Ratio), 
               alpha = 0.7) +
    geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
               size = 1, linetype = "longdash") + 
    ylab(ymalabel) + xlabel + 
    coord_cartesian(xlim = myXLimMA, ylim = myYLim) + 
    #colScale + 
    annotation_custom(tableGrob(rm_dat,theme=tt), 
                      #                                         gpar.corefill = gpar(fill = "grey85",
                      #                                                              col = "white"), 
                      #                                         gpar.rowfill = gpar(fill = "grey80",
                      #                                                             col = "white"),
                      #                                         gpar.colfill = gpar(fill = "grey80",
                      #                                                             col = "white")), 
                      #xmin = quantile(maData$A,probs=0.25),
                      #xmax = max(maData$A),
                      ymin = (myYLim[2]) - 0.25*myYLim[2], 
                      ymax = myYLim[2]) + 
    scale_y_continuous(breaks = seq(myYLim[1],myYLim[2],1))+ theme_bw()+
    theme( legend.justification = c(1,0),legend.position=c(1,0)) +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
print(maPlot)







