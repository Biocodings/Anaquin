#
# R script generated automatically for sequins.
#

x <- c(%1%)
y <- c(%2%)
z <- c(%3%)

m <- lm(y~x)

plot(x, y, xlab='Known Abundance', ylab='Measured Abundance')
