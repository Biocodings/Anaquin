#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

.transformPlot <- function(p)
{
    p <- p + theme(axis.title.x=element_text(face='bold', size=12))
    p <- p + theme(axis.title.y=element_text(face='bold', size=12))
    p <- p + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
    p <- p + theme(plot.title = element_text(size=15, face="bold"))
    p <- p + theme(legend.title = element_text(face="bold"))
    p <- p + theme(legend.key = element_blank())

    #
    # Calculating the aspect ratio to make the plot a square.
    #
    
    build <- ggplot_build(p)
    
    minX <- build$panel$ranges[[1]]$x.range[[1]]
    maxX <- build$panel$ranges[[1]]$x.range[[2]]
    minY <- build$panel$ranges[[1]]$y.range[[1]]
    maxY <- build$panel$ranges[[1]]$y.range[[2]]
    
    xrange <- maxX - minX
    yrange <- maxY - minY

    p <- p + coord_fixed(ratio=xrange/yrange)
 
    return (p)    
}