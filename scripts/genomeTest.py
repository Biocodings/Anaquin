#!/usr/bin/python

#
# This script provides implementation for comparing genomics regions.
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

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def countCoverage(bam, bed, binSize=9):
    
    r = bam.count_coverage(bed.chromID, bed.start, bed.end)
    
    x = []
    
    for i in range(len(r[0])):
        t = r[0][i] + r[1][i] + r[2][i] + r[3][i]
        x.append(t)

    bins = list(chunks(x, binSize))
    
    if (bins == 0):
        raise Exception('Failed to bin. No bin detected!')
    
    # The bins we're actually wanting to get
    y = []
    
    for i in range(len(bins)):
        t = []        
        for j in range(len(bins[i])):
            t.append(bins[i])
        y.append(numpy.mean(t))

    # Returns the coverage for each base
    return y

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

    print 'Testing: ' + name + ' ' + sign + ' ' + str(fold) + ' ' + str(stats)

    test.write(name + '\t' + str(numpy.mean(x)) + '\t' + str(numpy.std(x)) + '\t' + str(len(x)) + '\t' + str(stats) + '\t' + str(pval) + '\t' + str(fold) + '\t' + sign + '\n')

#
# Usage: <Alignment File.bam> <Region of Interest.bed> <Genome Background.bed>
#
#    Eg: python genomeTest.py NA12878_merged.Hg19.sorted.bam chrF.bed genome_wide_exon_regions.bed
#        python genomeTest.py filtered.bam custom.bed genome.bed
#        python genomeTest.py LLA072.unjoin_bwa.sort.bam unjoin_ACD_2.bed unjoin_B_2.bed
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

    logF = open('results.log', 'w')
    test = open('results.stats', 'w')    

    logF.write('Name\tLength\tCoverage\n')
    test.write('Name\tMean\tSD\tBases\tStats\tPValue\tFold\tSignificant\n')
    
    # Samples for the genome coverage
    x = []
    
    # Total number of bins for the genome
    xSize = 0
    
    binSize = 100
    
    for genome in genomes:

        # Length of the region
        size = genome.end - genome.start
        
        # Get the coverage in the region
        r = countCoverage(bam, genome, binSize)

        logF.write(genome.name + '\t' + str(size) + '\t' + str(sum(r)) + '\n')
        
        # Eg: [53.44444, 158.66666666, 264.111111,]
        x = x + r
        
        # Eg: 0 -> 9
        xSize = xSize + len(r)
        
    # Average genome average per bin...
    aver = sum(x) / xSize
    
    print ('\nNumber of bins is: ' + str(xSize))
    print ('Genome region is: ' + str(aver) + ' per bin' + '\n')

    #
    # Now, we'll apply statistical testing for each of the region we're interested
    #

    for custom in customs:
        
        # Get the coverage in the custom region
        r = countCoverage(bam, custom)
        
        # Let's compare it to the genomic average. Is the region statistically different?
        tTest(r, aver, custom.name, logF, test)

    print('\nGenerated results.log.')
    print('Generated results.stats.')
    
    bam.close()
    logF.close()
    test.close()    
