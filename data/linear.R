#
# R script generated automatically for sequin analysis.
#

x <- c(%1%)
y <- c(%2%)
z <- c(%3%)

# Color for each point
cols <- c(%4%)

# Pearson's correlation
p <- cor(x, y)

# Fit a simple linear regression model
m <- lm(y~x)

# Prints a summary statistics of the model
summary(m)

# Render the plot with filled circles
par(pch=19)

# Consturct a plot with some space on the top-corners for legends
plot(x, y, ylim=c(min(y), 2.0 * max(y)), col=cols, xlab='Log2 spike amount (attomoles/ul)', ylab='Log2 measured coverage (%5%)')

# Draw a line for LOS
abline(v=%6%, lty=2)

# Legend for correlation
legend("top", legend=c(paste('r:', round(cor(x,y), 4))), bty = "n")

# Legend for R2
legend("topright", legend=c(paste('R2:', round(summary(m)$r.squared, 4))), bty = "n")