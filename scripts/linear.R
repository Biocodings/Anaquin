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

# Render the plot with filled circles
par(pch=19)

# Consturct a plot with some space on the top-corners for legends
plot(seq.x, seq.y, xlim=c(min(x), 2.0 * max(x)), ylim=c(min(y), 2.0 * max(y)), col=seq.cols, xlab='Log2 spike amount (attomoles/ul)', ylab='Log2 measured coverage (%5%)')

# Draw a line for LOS
abline(v=%6%, lty=2)