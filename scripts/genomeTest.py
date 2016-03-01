#!/usr/bin/python

#
# This script provides implementation for comparing genomic regions.
#

import os
import sys
import pysam
import numpy
import scipy
import math

from statsmodels.stats.weightstats import ttest_ind

class BedLine:
    
    # Name of the chromosome
    chromID = None
    
    # Starting position
    start = -1
    
    # Ending position
    end = -1
    
    # Name of the region
    name = None

def readBed(file):
 
    #   
    # chrF  31539994    31540591    C_01_D
    #
    
    beds = []

    for line in open(file):
        line = line.strip()
        toks = line.split('\t')

        bed = BedLine()
        
        bed.name    = toks[3]
        bed.chromID = toks[0]
        bed.start   = int(toks[1])
        bed.end     = int(toks[2])
        
        beds.append(bed)
    
    return (beds)

def countCoverage(bam, bed):
    
    r = bam.count_coverage(bed.chromID, bed.start, bed.end)
    x = []
    
    for i in range(len(r[0])):
        t = r[0][i] + r[1][i] + r[2][i] + r[3][i]
        x.append(t)

    # Returns the coverage for each base
    return x

# Perform a one-sample t-test where x is the region of interest and n is the number of bases in the region
def tTest(x, aver, name, logF, test):

    t = scipy.stats.ttest_1samp(x, popmean=aver)
    
    stats = abs(t[0])
    pval  = t[1]
    
    if math.isinf(stats) or math.isinf(pval):
	print 'Warning: no coverage detected in: ' + name + ' .Skipped.'
        return

    #
    # <Region Name>  <Sample Mean>  <Sample Variance>  <Bases>  <Test Statistic>  <P-value>  <Fold>  <Signfiicant>
    #
    
    sign = 'F'
    
    if pval <= 0.05:
        sign = 'T'

    fold = numpy.mean(x) / aver

    print 'Testing: ' + name + ' ' + sign + ' ' + str(fold)

    test.write(name + '\t' + str(numpy.mean(x)) + '\t' + str(numpy.std(x)) + '\t' + str(len(x)) + '\t' + str(stats) + '\t' + str(pval) + '\t' + str(fold) + '\t' + sign + '\n')

#
# Usage: <Alignment File.bam> <Region of Interest.bed> <Genome Background.bed>
#
#    Eg: python genomeTest.py NA12878_merged.Hg19.sorted.bam chrF.bed genome_wide_exon_regions.bed
#        python genomeTest.py filtered.bam custom.bed genome.bed
#
#          Generates:
#
#               - results.log
#               - results.stats (statistical results)
#

if __name__ == '__main__':

    if (len(sys.argv) != 4):
        print 'GenomeTest is a program for performing statistical test of a chosen region to the genome. The user will need to choose the genome'
        print 'and regions of interest.\n\nThe script performs a one-sample t-test, one tailed for each of the region of interests against the genome.\n'
        print '     Usage: <BAM file> <Region of interest.bed> <Genome background.bed>'
        print '     Outputs: results.log and results.stats.\n\nresults.log shows the detailed logging whereas results.stats shows the statistical tests\n'
        print 'Eg: genomeTest.py /usr/bin/genomeTest/filtered.bam /usr/bin/genomeTest/custom.bed /usr/bin/genomeTest/genome.bed\n'
        sys.exit()

    align  = sys.argv[1]
    inter  = sys.argv[2]
    genome = sys.argv[3]

    print ('Region interest:   ' + inter)
    print ('Genome background: ' + genome)

    # Region of interest
    customs = readBed(inter)
    
    # Genomic region
    genomes = readBed(genome)

    bam = pysam.AlignmentFile(align, 'rb')

    print 'Calculating the average for the genomic regions...'

    #
    # Calculate the average coverage for the genomic region
    #

    n = 0
    total = 0
    
    logF = open('results.log', 'w')
    test = open('results.stats', 'w')    

    logF.write('Name\tLength\tCoverage\n')
    test.write('Name\tMean\tSD\tBases\tStats\tPValue\tFold\tSignificant\n')
    
    for genome in genomes:

        # Length of the region
        size = genome.end - genome.start
        
        # Get the coverage in the region
        r =  countCoverage(bam, genome)

        logF.write(genome.name + '\t' + str(size) + '\t' + str(sum(r)) + '\n')
        
        n = n + size
        total = total + sum(r)

    aver = total / n

    print ('\nGenome region is: ' + str(aver) + ' per base (average)' + '\n')

    #
    # Now, we'll apply statistical testing for each of the region we're interested
    #

    for custom in customs:
        
        # Get the coverage in the custom region
        x = countCoverage(bam, custom)
        
        # Let's compare it to the genomic average. Is the region statistically different?
        tTest(x, aver, custom.name, logF, test)

    print('Generated results.log. Please check it for logging.')
    print('Generated results.stats. Please check it for statistics.')
    
    bam.close()
    logF.close()
    test.close()    
