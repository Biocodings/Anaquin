#
# R script generated automatically for sequin analysis.
#

x <- c(%1%)
y <- c(%2%)
z <- c(%3%)

# Pearson's correlation
p <- cor(x, y)

# Fit a simple linear regression model
m <- lm(y~x)

# Prints a summary statistics of the model
summary(m)

plot(x, y, xlab='Log Known Abundance (attomoles/ul)', ylab='Log Measured Abundance (k-mer average)')