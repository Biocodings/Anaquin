#
# Anaquin - Sequin statistical analysis. Version 1.0.
#
# This R script was generated at %1%.
#
#    %2%
#

library('ggplot2')

# Convert a linear model to string
lm_eqn <- function(d)
{
    m <- lm(y ~ x, d);
    eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

d <- data.frame(x=c(%3%),
                y=c(%4%),
                ids=c(%5%))

p <- ggplot(data = d, aes(x = x, y = y))
p <- p + xlab('Log2 spike amount (attomoles/ul)')
p <- p + ylab('Log2 measured coverage (%6%)')
p <- p + geom_point()
p <- p + geom_smooth(method = "lm", se=FALSE, formula = y ~ x)
p <- p + geom_text(x = min(d$x), y = max(d$y), label = lm_eqn(d), parse = TRUE, hjust=0)
print(p)

# Pearson's correlation
#p <- cor(x, y)

# Fit a simple linear regression model
#m <- lm(y~x)

# Prints a summary statistics of the model
#summary(m)