#!/usr/bin/python

#
# This script constructs a distribution plot for coverage per base for hg38 and chrT.
#

import os
import sys

def run(cmd):
    os.system(cmd)

if __name__ == '__main__':

    file = sys.argv[1]


    run('samtools view -bS aligned.sam > aligned.bam')
    run('samtools sort aligned.bam aligned_orderd.bam')
    run('bedtools genomecov -bg -ibam aligned_ordered.bam > aligned.bedgraph')

    run('grep chrT aligned.bedgraph > chrT.bedgraph')
    run('grep -v chrT aligned.bedgraph > hg38.bedgraph')

    
    
