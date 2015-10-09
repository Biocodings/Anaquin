#!/usr/bin/python

#
# This script constructs distribution for base coverage for chrT and hg38.
#

import os
import sys

def run(cmd):
    os.system(cmd)

if __name__ == '__main__':

    # The alignment file
    file = sys.argv[1]

    
    
    
    
    # run('bedtools genomecov -bg -ibam aligned_chrT.bam > aligned_chrT.bedgraph')
    # run('bedtools genomecov -bg -ibam aligned_sampled.bam > aligned_sampled.bedgraph')
    