## Spike - Garvan Institute
##
## Copyright (C) 2015   Dr Timothy Mercer & Ted Wong
##
##  RQuantLib is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 2 of the License, or
##  (at your option) any later version.
##
##  RQuantLib is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with Spike. If not, see <http://www.gnu.org/licenses/>.



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










