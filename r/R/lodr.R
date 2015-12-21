#
#  Copyright (C) 2015 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#


m <- read.csv('/Users/tedwong/Desktop/LODR_data_TED.csv', row.names=1)
m$expected.log2FoldChange = as.factor(abs(m$expected.log2FoldChange))
m <- m[!is.na(m$padj),]
plotLODR(m$baseMean, m$padj, m$expected.log2FoldChange)

plotLODR <- function(counts, pvals, ratios, cutoff=0.01, xname='Average Counts', yname='DE Test P-values')
{
    xrange <- c(1, 28199)
    
    d <- data.frame(x=counts, y=pvals, ratios=ratios)
    
    ggplot(d, aes(x=x, y=y,colour=ratios)) + 
           geom_point(size = 6) +
           xlab(xname) +
           ylab(yname) + 
           
           scale_x_log10(limits = xrange) + 
           scale_y_log10(breaks = c(1e-300,1e-200,1e-100,1e-10, 1.00)) +
        
           #geom_ribbon(data = lineDat, aes(x = x.new, 
           #                                y = fitLine, 
           #                                ymin=fitLower, 
           #                                ymax=fitUpper,
           #            fill = Ratio), alpha = 0.3,
           #colour=NA,show_guide=FALSE) + 
        
           #geom_line(data = lineDat,aes(x = x.new, y=fitLine, 
           #          colour = Ratio),show_guide = FALSE) +

           #colScale + fillScale +

           geom_hline(yintercept = cutoff, linetype = 2, size = 2 ) + 
    
           theme_bw()
    
           #geom_segment(data = arrowDat, aes(x = x,
           #                                     y = y, 
           #                                  xend = xend,
           #                                     yend = yend, 
           #                                  colour = Ratio), 
           #             lineend = "round", arrow=arrow(length=unit(0.5,"cm")), 
           #             size = 2, alpha = 0.6) +
}
