#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#


library(grid)

arrowDat <- read.csv('arrowDat', row.names=1)
lineDat <- read.csv('lineDat', row.names=1)
pval.res <- read.csv('pval.res', row.names=1)
xrange <- c(1, 28199)
pval.cutoff <- 0.00679406


LODRplot <- ggplot(pval.res, aes(x=MnSignal,
                                 y=Pval,colour=Ratio)) + 
    geom_point(size = 6) + 
    scale_x_log10(limits = xrange) + 
    scale_y_log10(breaks = c(1e-12,1e-10,1e-8,1e-5,1e-4,
                             1e-3,1e-2,1e-1,1e0))+
    geom_ribbon(data = lineDat, aes(x = x.new, 
                                    y = fitLine, 
                                    ymin=fitLower, 
                                    ymax=fitUpper,
                                    fill = Ratio), alpha = 0.3,
                colour=NA,show_guide=FALSE) + 
    geom_line(data = lineDat,aes(x = x.new,
                                 y=fitLine, 
                                 colour = Ratio),show_guide = FALSE) + 
    #colScale + fillScale + xlabDE + ylab("DE Test P-values") + 
    geom_hline(yintercept = pval.cutoff, linetype = 2, size = 2 ) + 
    geom_segment(data = arrowDat, aes(x = x,
                                      y = y, 
                                      xend = xend,
                                      yend = yend, 
                                      colour = Ratio), 
                 lineend = "round", arrow=arrow(length=unit(0.5,"cm")), 
                 size = 2, alpha = 0.6) + theme_bw() 
#+ legendPos  

LODRplot

