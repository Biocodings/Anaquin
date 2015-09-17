#!/usr/bin/python

#
# This script constructs a distribution plot for coverage per base for hg38 and chrT.
#

import os
import sys

def run(cmd):
    print('---------------------------------------------------------------')
    print cmd
    print('---------------------------------------------------------------')
    #os.system(cmd)

if __name__ == '__main__':

    file = sys.argv[1]

    # This is needed to estimate the subsampling ratio
    #run('grep -v chrT ' + file + '| wc')

    # Number of reads in the alignment excluding chrT
    reads = 150586265 # From the wc command above

    # This is the fraction we'll need to subsample the genome
    ratio = 26011845 / reads

    ratio = 0.1727372    
    
    #
    # Eg: grep '@\|chrT' aligned.sam > temp.sam
    #
    run('grep \'@\|chrT\' ' + file + ' > temp.sam')

    #
    # Eg: grep -v chrT aligned.sam | grep -v '*'
    #
    run('grep -v chrT aligned.sam | grep -v \'*\'')

    #
    # Eg: grep -v chrT aligned.sam | grep -v '*' | samtools view -S -s 0.05 - >> temp.sam
    #
    run('grep -v chrT aligned.sam | grep -v '*' | samtools view -S -s 0.05 - >> temp.sam')
    
    #
    # Eg: samtools view -Sb  temp.sam | samtools sort - sorted
    #
    
    #
    # Eg: bedtools genomecov -bg -ibam sorted.bam > sorted.bedgraph
    #

    #
    # Eg: grep -v chrT sorted.bedgraph > hg38.bedgraph
    # Eg: grep chrT sorted.bedgraph > chrT.bedgraph
    #







    # run('grep -v chrT aligned.sam > hg38.sam')
    #   run('grep chrT aligned.sam > chrT.sam')
    #  run('samtools view -bS chrT.sam > chrT.bam')
    
    # run('wc hg38.sam') # How many reads? 150586265
    #run('wc chrT.sam') # How many reads? 150586265
    
    #
    # It's probably too time-consuming to construct a distribution for the genome. Why not just take a sufficent random samples?
    #
    
    
    
    #run('samtools view -s 0.1727372 -b chrT.bam > sampled.bam')

    #    run('samtools sort chrT.bam aligned_chrT')
    #   run('samtools sort sampled.bam aligned_sampled')
    #  run('bedtools genomecov -bg -ibam aligned_chrT.bam > aligned_chrT.bedgraph')
    # run('bedtools genomecov -bg -ibam aligned_sampled.bam > aligned_sampled.bedgraph')

    #
    # Draw each individual bedgraph on the same plot...
    #

    #    run('grep chrT aligned.bedgraph > chrT.bedgraph')
    #    run('grep -v chrT aligned.bedgraph > hg38.bedgraph')