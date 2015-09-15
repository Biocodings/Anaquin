#!/usr/bin/python

#
# This script constructs a distribution plot for coverage per base for hg38 and chrT.
#

import os
import sys

def run(cmd):
    os.system(cmd)

if __name__ == '__main__':

    run('grep chrT aligned.sam > chrT.sam')
    run('grep -v chrT aligned.sam > hg38.sam')

    run('samtools view -bS chrT.sam > chrT.bam')
    run('samtools view -bS hg38.sam > hg38.bam')

    run('samtools sort chrT.bam chrT_sorted')
    run('samtools sort hg38.bam hg38_sorted')    
    run('rm chrT.bam')
    run('rm hg38.bam')        

    # Get info.csv from Dropbox...

    run('bedtools genomecov -bg -ibam chrT_sorted.bam > chrT_sorted.bedgraph')
    run('bedtools genomecov -bg -ibam hg38_sorted.bam > hg38_sorted.bedgraph')

    # How to draw a distribution plot from bedgraph??
    
    
    
