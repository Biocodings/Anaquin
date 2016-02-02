#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at 03-02-2016 12:50:47.
#
#    ./anaquin -t VarDiscover -rbed data/VARQuin/AVA017.v032.bed -rvcf data/VARQuin/AVA009.v032.vcf -m data/VARQuin/MVA012.v013.csv -soft VarScan -ufiles varscan.tab 
#

library(Anaquin)

# Change this for your significance level
signLevel <- 0.10

# Read the data file for true positives
tps <- read.csv('/Users/tedwong/Sources/QA/output/VarDiscover_TP.csv', sep='\t')

# Read the data file for false positives
fps <- read.csv('/Users/tedwong/Sources/QA/output/VarDiscover_FP.csv', sep='\t')

# Filter only the significant FPs
fps <- fps[fps$PValue < signLevel,]

tps$label <- 'TP'
fps$label <- 'FP'

data <- rbind(tps, fps)

data$name <- paste(data$Sequin, data$Position, sep='_')
data$name <- paste(data$name, data$Type, sep='_')

# Create a TransQuin data set for Anaquin
data <- VarQuin(seqs=data$name, pval=data$PValue, rRead=data$RefRead, vRead=data$VarRead, eAFreq=data$EAlleleF, label=data$label, type=data$Type)

plotROC.VarQuin(data)
