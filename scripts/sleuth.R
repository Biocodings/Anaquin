#
# This script runs differential analysis with k-mers with the sleuth R-package.
#

library("sleuth")

# Where the kallisto results are
abunds <- list(%1%)

# Where the kallisto paths are
paths <- sapply(abunds, dirname)

s2c <- data.frame(sample=c(%2%), condition=c(%3%))
s2c <- dplyr::mutate(s2c, path=paths)

so <- sleuth_prep(s2c, ~condition)
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'conditionB')

results <- sleuth_results(so, 'conditionB')

# This is the CSV file we'll need to parse
write.csv(results, file='%4%', row.names=FALSE, quote=FALSE)
 