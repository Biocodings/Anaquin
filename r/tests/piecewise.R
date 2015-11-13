
d <- read.csv('/Users/tedwong/Desktop/RXMA_guided_obsExp.csv', row.names=1)
colnames(d) <- c('x', 'y')
d <- d[d$x<10000,]
d <- d[d$y>0,]

plot(log2(d$x), log2(d$y))

r <- softLimit(log2(d$x), log2(d$y))

