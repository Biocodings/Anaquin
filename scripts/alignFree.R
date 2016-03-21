
# Mapping for transformed k-mers
kmerMap <- '/Users/tedwong/Desktop/kTest/index.mapping'

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

r <- aggregate(kmerSeq$abund, list(kmerSeq$seq), FUN=mean)

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
merged$measred <- merged$m_abund1 / merged$m_abund2




