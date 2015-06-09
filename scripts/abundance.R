#
# R script generated automatically for sequins.
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

plot(x, y, xlab='Known Abundance', ylab='Measured Abundance')
