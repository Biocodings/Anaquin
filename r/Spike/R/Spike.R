# Spike - Garvan institute 
#
#



RNASequins <- function()
{
    ms <- list("A" = ("/Users/tedwong/Sources/QA/data/RNA/mixture_A.csv"),
               "B" = ("/Users/tedwong/Sources/QA/data/RNA/mixture_B.csv"))
    ms
}

SP_DESeq <- function(object, test = c("Wald", "LRT"), fitType = c("parametric", 
                                                                  "local", "mean"), betaPrior, full = design(object), reduced, 
                     quiet = FALSE, minReplicatesForReplace = 7, modelMatrixType, 
                     parallel = FALSE, BPPARAM = bpparam())
{
    # We can do the origianl test... But what to show?
    assay(object)[1:3,]
}




#library("airway")
#data("airway")
#se <- airway
#library("DESeq2")
#ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
#d <- DESeq(ddsSE)










