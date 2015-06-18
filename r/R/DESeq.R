#  Copyright (C) 2015 - Garvan Institute (Dr Timothey Mercer, Dr Wendy Chen, Ted Wong)
#
#  Anaquin is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  Anaquin is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Anaquin If not, see <http://www.gnu.org/licenses/>.

#
# This script demostrates a common workflow of DESeq2 and Anaquin
#

library('DESeq2')

# Load the column data that describes the experiment
col <- read.csv('/Users/tedwong/Sources/QA/r/data/colData.csv')

# Load the count matrix
data <- read.csv('/Users/tedwong/Sources/QA/r/data/airway.csv', row.names=1)

# Construct an object for DESeq2
se <- DESeqDataSetFromMatrix(as.matrix(data), colData=col, design=~dex)

# Run the differential expression test
dds <- DESeq(se)

# Gather the results
res <- results(dds)

# Call to Anaquin
RNAAnalyze(res)




