library(data.table)

# Mapping for transformed k-mers
kmerMap <- '/Users/tedwong/Sources/QA/data/FusQuin/AFU011.v032.map'

# Generated abundance
kmerAbund <- '/Users/tedwong/Desktop/kTest/output/matched.kmers'

# Information about the k-mers in the sequins
kmerSeq <- '/Users/tedwong/Sources/QA/data/FusQuin/AFU010.v032.csv'


kmerMap <- read.table(kmerMap, header=FALSE)
colnames(kmerMap) <- c('orig', 'trans')

kmerAbund <- read.table(kmerAbund, header=FALSE)
colnames(kmerAbund) <- c('name', 'abund')

kmerSeq <- read.table(kmerSeq, header=FALSE)
colnames(kmerSeq) <- c('seq', 'id', 'status', 'start', 'end', 'kmer')
kmerSeq <- kmerSeq[c(1,2,3,6)]

#
# 1: Fitler only the k-mers spanning the fusion points (the file would have everything)
#

kmerSeq <- kmerSeq[kmerSeq$status=='FS' | kmerSeq$status=='NS',]

#
# 2: Define a function that we can use to map k-mers
#

orig2trans <- function(kmer)
{
    return (kmerMap[kmerMap$orig %in% kmer,]$trans)
}

trans2orig <- function(kmer)
{
    return (kmerMap[kmerMap$trans %in% kmer,]$orig)
}

#
# 3: Assign abundance to the kmerSeq table. 
#

f1 <- function(x)
{
    return (kmerAbund[kmerAbund$name %in% orig2trans(x["kmer"]),]$abund)
}

kmerSeq$abund <- apply(kmerSeq, 1, f1)

#
# 4: Aggregate into each individual sequin
#

r <- aggregate(kmerSeq$abund, list(kmerSeq$seq), FUN=sum)

#
# 5: Read mixture file for sequins
#

mix <- loadMixture.FusQuin()
merged <- normalFusion(mix)

merged$m_abund1 <- NA
merged$m_abund2 <- NA

for (i in 1:nrow(r))
{
    row <- r[i,]
    a <- row$Group.1

    if (nrow(merged[merged$names1==a,]) > 0)
    {
        merged[merged$names1==a,]$m_abund1 <- row$x
    }
    else if (nrow(merged[merged$names2==a,]) > 0)
    {
        merged[merged$names2==a,]$m_abund2 <- row$x
    }
    else
    {
        stop ('????')
    }
}

merged <- merged[!is.na(merged$m_abund1),]
merged$measured <- merged$m_abund1 / merged$m_abund2


#> merged[log2(merged$measured) <= -5,]
#names1    names2    abund1   abund2     abund m_abund1 m_abund2      measred     measured
#4  NG1_10_P2 FG1_10_P2  483.4005   7.5570 63.967249        7    73130 9.571995e-05 9.571995e-05
#10  NG1_2_P2  FG1_2_P2 1933.5938 483.4005  3.999983       45     9268 4.855416e-03 4.855416e-03

a <- merged[log2(merged$measured) >= -5,]
x <- log2(a$abund)
y <- log2(a$measured)


#kmerSeq[kmerSeq$seq=='NG1_10_P2'|kmerSeq$seq=='FG1_10_P2',]

#library(ggplot2)

#p <- ggplot(data=a, aes(x=m_abund1, y=m_abund2)) +
#    geom_point() +
#    geom_smooth(method='lm', formula=y~x)            +
#    theme_bw()
#print(p)





