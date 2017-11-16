
leftStargazer <- function(x)
{
    collapse <- function(st) paste(st, collapse="")
    st <- gsub(collapse(rep("c", 3)), collapse(rep("l",3)), x)
    cat(st[4:length(st)])
}

#
# Reference files
#

refs.d <- c('Reference sequin variants', 'Sequin regions (hg38)')
refs.f <- c('/path/to/sequins_variants.hg38.vcf', '/path/to/L.V.97.human.raw_variants.vcf')
refs.s <- c('Default', 'Default')
refs   <- data.frame('Description'=refs.d, 'File Name'=refs.f, 'Source'=refs.s)

#
# Input files
#

inputs.d <- c('User sequin variants', 'User sample variants', 'User regions')
inputs.f <- c('L.V.97.sequin.raw_variants.vcf', 'L.V.97.human.raw_variants.vcf', 'regions.bed')
inputs.t <- c('VCF annotation', 'VCF annotation', 'BED annoation')
inputs   <- data.frame('Description'=inputs.d, 'File Name'=inputs.f, 'Type'=inputs.t)

#
# Germline variants
#

germs.f <- c('All variants', 'Comparison regions', 'Comparison size (bp)')
germs.a <- c('CS,CI,GS,GI,LoGC,VLGC,HiGC,VHGC,MS,SH,LH,SD,LD,ST,LT,SQ,LQ', '270', '157134')
germs   <- data.frame('Description'=germs.f, 'Information'=germs.a)

#
# Germline performance
#

germs.c <- c('Total', 'SNVs', 'Indels')
germs.r <- c('Reference sequin variants', 'True positives (TP)', 'False positives (FP)', 'False negatives (FN)', 'Sensitivity', 'Precision', 'F1 Score', 'FDR Rate', 'FPR (FP/kb)')
germs.t <- c(170, 152, 21, 21, 0.89, 0.873, 0.883, 0.126, 0.13)
germs.s <- c(170, 152, 21, 21, 0.89, 0.873, 0.883, 0.126, 0.13)
germs.i <- c(170, 152, 21, 21, 0.89, 0.873, 0.883, 0.126, 0.13)
germs.m <- data.frame(germs.t, germs.s, germs.i)
colnames(germs.m) <- germs.c
rownames(germs.m) <- germs.r
