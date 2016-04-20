#
#  Copyright (C) 2016 - Garvan Institute of Medical Research
#
#  Ted Wong, Bioinformatic Software Engineer at Garvan Institute
#

library('RUnit')
library('Anaquin')

test_3 <- function()
{
    #
    # Experiment with three replicates at the exon level
    #
    
    data <- read.csv('/Users/tedwong/Desktop/sequins_exons_full_TED.csv', row.names=1)
    colnames(data) <- c('X', 'A1', 'A2', 'A3')
    
    data$X  <- as.numeric(as.character(data$X))
    data$A1 <- as.numeric(as.character(data$A1))
    data$A2 <- as.numeric(as.character(data$A2))
    data$A3 <- as.numeric(as.character(data$A3))
    
    data$Mean <- rowMeans(data[,c(-1)])
    head(data)
    
    A1 <- data.frame(X=data$X, Y=data$A1)
    A2 <- data.frame(X=data$X, Y=data$A2)
    A3 <- data.frame(X=data$X, Y=data$A3)
    r  <- data.frame(X=c(A1$X, A2$X, A3$X), Y=c(A1$Y, A2$Y, A3$Y))
    head(r)
    
    # Option 1: X vs Mean?
    plot(log10(data$X), log10(data$Mean))
    
    # Option 2: X vs All?
    plot(log10(r$X), log10(r$Y))
    
    r <- r[r$X!=0,]    
    r <- r[r$Y!=0,]
    
    # Let's use option 2... More data points...
    r <- plotInflection(log10(r$X), log10(r$Y), showDetails=TRUE)
    
    # Breakpoints and R2s
    r$breaks
    
    # Regression for the left
    summary(r$model$lModel)
    
    # Regression for the right
    summary(r$model$rModel)
}

test_2 <- function()
{
    #
    # Experiment with three replicates at the isoform level
    #
    
    data <- read.csv('/Users/tedwong/Desktop/sequins_isoforms_full_TED.csv', row.names=1)
    data <- data[-c(1),]
    colnames(data) <- c('X', 'A1', 'A2', 'A3')
    
    data$X  <- as.numeric(as.character(data$X))
    data$A1 <- as.numeric(as.character(data$A1))
    data$A2 <- as.numeric(as.character(data$A2))
    data$A3 <- as.numeric(as.character(data$A3))
    
    data$Mean <- rowMeans(data[,c(-1)])
    head(data)
    
    A1 <- data.frame(X=data$X, Y=data$A1)
    A2 <- data.frame(X=data$X, Y=data$A2)
    A3 <- data.frame(X=data$X, Y=data$A3)
    r  <- data.frame(X=c(A1$X, A2$X, A3$X), Y=c(A1$Y, A2$Y, A3$Y))
    head(r)
    
    # Option 1: X vs Mean?
    plot(log10(data$X), log10(data$Mean))
    
    # Option 2: X vs All?
    plot(log10(r$X), log10(r$Y))
    
    r <- r[r$X!=0,]    
    r <- r[r$Y!=0,]
    
    summary(lm(log10(r$Y) ~ log10(r$X)))
    
    # Let's use option 2... More data points...
    r <- plotInflection(log10(r$X), log10(r$Y), showDetails=TRUE)
    
    # Breakpoints and R2s
    r$breaks
    
    # Regression for the left
    summary(r$model$lModel)
    
    # Regression for the right
    summary(r$model$rModel)
}

test_1 <- function()
{
    #
    # Experiment with three replicates at the gene level
    #
    
    data <- read.csv('/Users/tedwong/Desktop/sequins_genes_full_TED.csv', row.names=1)
    data <- data[-c(1),]
    colnames(data) <- c('X', 'A1', 'A2', 'A3')
    
    data$X  <- as.numeric(as.character(data$X))
    data$A1 <- as.numeric(as.character(data$A1))
    data$A2 <- as.numeric(as.character(data$A2))
    data$A3 <- as.numeric(as.character(data$A3))

    data$Mean <- rowMeans(data[,c(-1)])
    head(data)
    
    A1 <- data.frame(X=data$X, Y=data$A1)
    A2 <- data.frame(X=data$X, Y=data$A2)
    A3 <- data.frame(X=data$X, Y=data$A3)
    r  <- data.frame(X=c(A1$X, A2$X, A3$X), Y=c(A1$Y, A2$Y, A3$Y))
    head(r)

    # Option 1: X vs Mean?
    plot(log10(data$X), log10(data$Mean))

    # Option 2: X vs All?
    plot(log10(r$X), log10(r$Y))

    summary(lm(log10(r$Y) ~ log10(r$X)))

    # Let's use option 2... More data points...
    r <- plotInflection(log10(r$X), log10(r$Y), showDetails=TRUE)
    
    # Breakpoints and R2s
    r$breaks

    # Regression for the left
    summary(r$model$lModel)

    # Regression for the right
    summary(r$model$rModel)
}

test_4 <- function()
{
    #
    # Differential analysis: expected abundance vs measured abundance
    #
    
    data <- TransQuin(seqs     = c('R1_101','R1_102','R1_103','R1_11','R1_12','R1_14','R1_22','R1_23','R1_24','R1_31','R1_32','R1_33','R1_41','R1_43','R1_51','R1_52','R1_53','R1_61','R1_62','R1_72','R1_73','R1_81','R1_82','R1_83','R1_91','R1_92','R1_93','R2_1','R2_105','R2_115','R2_116','R2_117','R2_150','R2_151','R2_152','R2_153','R2_154','R2_20','R2_24','R2_27','R2_28','R2_32','R2_33','R2_37','R2_38','R2_41','R2_42','R2_45','R2_46','R2_47','R2_53','R2_54','R2_57','R2_59','R2_6','R2_60','R2_63','R2_65','R2_67','R2_68','R2_7','R2_71','R2_72','R2_73','R2_76'),
                      expected = c(0.125000, 0.062500, 0.500000, 0.062500, 2.000000, 0.500000, 0.062500, 1.000000, 0.062500, 8.000000, 0.062500, 0.500000, 0.500000, 16.000000, 0.062500, 16.000000, 2.000000, 0.125000, 0.250000, 0.062500, 0.125000, 1.000000, 1.000000, 4.000000, 16.000000, 0.500000, 8.000000, 16.000000, 8.000000, 4.000000, 1.000000, 1.000000, 0.125000, 1.000000, 1.000000, 0.500000, 0.500000, 1.000000, 1.000000, 0.125000, 4.000000, 16.000000, 2.000000, 0.250000, 4.000000, 0.250000, 2.000000, 1.000000, 2.000000, 8.000000, 1.000000, 1.000000, 0.125000, 1.000000, 8.000000, 16.000000, 1.000000, 1.000000, 0.062500, 0.062500, 1.000000, 1.000000, 0.250000, 2.000000, 8.000000),
                      measured = c(0.408303, 0.259021, 0.813703, 0.089496, 0.832248, 0.753653, 0.131169, 1.447562, 0.083715, 12.428689, 0.204704, 0.878076, 0.829980, 10.837457, 0.101030, 4.035239, 3.107969, 0.636319, 1.018794, 0.808029, 0.167796, 1.441933, 1.876315, 6.543644, 2.993838, 1.009377, 9.166627, 3.135811, 3.867229, 7.964597, 1.060898, 1.717270, 0.184347, 1.091200, 1.523789, 1.070212, 0.760873, 1.426753, 1.458052, 0.474776, 1.430541, 2.191430, 0.690605, 0.942656, 1.074324, 0.377179, 1.028790, 1.024330, 1.116275, 14.444278, 1.063578, 1.631630, 0.429272, 1.039425, 9.363022, 11.534607, 1.377474, 0.611819, 0.397484, 12.387706, 1.542621, 1.150802, 1.151967, 1.681149, 1.160848))
    r <- plotScatter(data, xname = 'Expected log2 fold change', yname = 'Measured log2 fold change')
    
    checkEquals(r$xname, 'Expected log2 fold change')
    checkEquals(r$yname, 'Measured log2 fold change')    
}

#test_1()
#test_4()