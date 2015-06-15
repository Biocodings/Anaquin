#
# R script generated automatically for sequin analysis.
#

seq.x <- c(%1%)
seq.y <- c(%2%)
seq.z <- c(%3%)

# Color for each point
seq.cols <- c(%4%)

# Pearson's correlation
p <- cor(x, y)

# Fit a simple linear regression model
m <- lm(y~x)

# Prints a summary statistics of the model
summary(m)

plot(seq.x, seq.y, col=seq.cols, xlab='Log Known Abundance (attomoles/ul)', ylab='Log Measured Abundance (k-mer average)')