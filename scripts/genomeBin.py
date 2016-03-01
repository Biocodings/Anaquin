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

from Bio import SeqIO

def readBAM(bam, chromID, start, end):
    
    try:
        r = bam.count_coverage(chromID, start, end)
    except:
        print 'Invalid ' + chromoID
        return None

    # Sum of coverage acroess the region
    return (sum(r[0]) + sum(r[1]) + sum(r[2]) + sum(r[3]))  / (end - start)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

#
# Usage: <Alignment File.bam> <Region of Interest.bed> <Genome Background.bed>
#
    #    Eg: genomeBin.py filtered.bam chrT.fa
#

if __name__ == '__main__':

    # The file for the genome (FASTA)
    genome = sys.argv[1]

    dict = {}

    fasta_sequences = SeqIO.parse(open(sys.argv[2]), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        dict[name] = sequence
        print 'Found: ' + name

    # Read in for the BAM file
    bam = pysam.AlignmentFile(sys.argv[1], 'rb')
    
    binSize = 600

    for chrID in dict:
        print 'Analyzing: ' + chrID
        
        # Split the chromosome into bins with equal size
        #bins = list(chunks(dict[chrID], binSize))        
        #print 'Number of bins: ' + str(len(bins))
        
        i = 0
        j = len(dict[chrID])

        # Frequency table
        freq = {}

        while i < j:

            # Read the coverage for bin
            r = readBAM(bam, chrID, i, i + binSize)

            if r is None:
                continue

            if r not in freq:
                freq[r] = 0            

            freq[r] = freq[r] + 1

            # Prepare for the next bin...
            i = i + binSize
        
        writer = open('results.csv', 'w')

        for i in freq:
            writer.write(chrID + '\t' + str(i) + '\t' + str(freq[i]) + '\n')

        writer.close()
