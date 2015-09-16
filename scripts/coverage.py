#!/usr/bin/python

#
# This script constructs a distribution plot for coverage per base for hg38 and chrT.
#

import os
import sys

def run(cmd):
    print cmd
    os.system(cmd)

if __name__ == '__main__':

    #file = sys.argv[1]

    run('grep -v chrT aligned.sam > hg38.sam')
    run('grep chrT aligned.sam > chrT.sam')
    run('samtools view -bS chrT.sam > chrT.bam')
    
    #
    # It's probably too time-consuming to construct a distribution for the genome. Why not just take a sufficent random samples?
    #
    
    # Number of reads in the alignment excluding chrT
    reads = 1000
    
    # This is the fraction we'll need to subsample the genome
    ratio = 0.05
    
    run('samtools view -s 0.25 -b chrT.bam > sampled.sam')
    run('samtools view -bS sampled.sam > sampled.bam')    

    run('samtools sort chrT.bam aligned_chrT.bam')
    run('samtools sort sampled.bam aligned_sampled.bam)
    
    run('bedtools genomecov -bg -ibam aligned_chrT.bam > aligned_chrT.bedgraph')
    run('bedtools genomecov -bg -ibam aligned_sampled.bam > aligned_sampled.bedgraph')

#    run('grep chrT aligned.bedgraph > chrT.bedgraph')
#    run('grep -v chrT aligned.bedgraph > hg38.bedgraph')

    
    
