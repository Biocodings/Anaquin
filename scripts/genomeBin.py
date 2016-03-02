#!/usr/bin/python

#
# This script provides implementation for genomical regions.
#

import os
import sys

class BedLine:
    
    # Name of the chromosome
    chromID = None
    
    # Starting position
    start = None
    
    # Ending position
    end = None
    
    # Name of the region
    name = None

class GenomeCoverage:

    # Name of the chromosome
    chromID = None

    # The position of the coverage
    start = None
    
    # Number of reads aligning to the position
    reads = None

def readGenome(file):
 
    gens = []

    for line in open(file):
        line = line.strip()
        
        # Eg: chr21  28  0
        toks = line.split('\t')

        gen = GenomeCoverage()
        
        gen.chromID = toks[0]
        gen.start   = int(toks[1])
        gen.reads   = int(toks[2])
    
        gens.append(gen)
    
    return gens

def readBed(file):
 
    beds = []

    for line in open(file):
        line = line.strip()
        
        # Eg: chrF  31539994    31540591    C_01_D
        toks = line.split('\t')

        bed = BedLine()
        
        bed.name    = toks[3]
        bed.chromID = toks[0]
        bed.start   = int(toks[1])
        bed.end     = int(toks[2])
        
        beds.append(bed)
        
    return beds

if __name__ == '__main__':

    mode = sys.argv[1]

    # This tool is used to filter per-base coverage results from bedtools by given region in a Bed file
    if mode == 'filterByBed':        
        bedFile = sys.argv[2]
        genFile = sys.argv[3]
        
        #print ('Reading the Bed file...')
        bed = readBed(bedFile)
        
        #print ('Reading the genome...')
        gen = readGenome(genFile)

        #print ('Number of genomes: ' + str(len(gen)))
        #print ('Number of beds: ' + str(len(bed)))

        for i in gen:
            if i.reads >= 5:
                for j in bed:
                    if (i.start >= j.start and i.start <= j.end):
                        #print str(i.start) + ' ' + str(j.start) + ' ' + str(j.end)
                        print (i.chromID + '\t' + str(i.start) + '\t' + str(i.reads))
                        break
        
    elif mode == 'countBedTools':
        file = sys.argv[2]        

        # Size of each bin
        binSize = 1

        i = 0
        n = 0

        min = sys.maxint
        max = -1
        dict = {}

        for line in open(file):
            line = line.strip()

            # Eg:  C_01_A  21  0
            toks = line.split('\t')

            # Coverage per base
            cov = float(toks[2])

            n = n + cov            
            i = i + 1
        
            if i == binSize:

                t = n
                n = int(n)
                
                if n not in dict:
                    dict[n] = 0
                    
                # We've just found this coverage, increment the counter by 1
                dict[n] = dict[n] + 1
                
                if n < min:
                    min = n
                if n > max:
                    max = n
                    
                i = 0
                n = 0

        #print 'Min: ' + str(min)
        #print 'Max: ' + str(max)
        #print 'Filling out the gaps...'
                
        for i in range(max):
            if i not in dict:
                dict[i] = 0

        # Print out the freqency table
        for key in sorted(dict):
            print (str(key) + '\t' + str(dict[key]))

    else:
        # What's the proportion of reads for the standards?
        ratio = 0.05
    
        # samtools idxstats NA12878_merged.Hg19.sorted.REF_chr21.bam
        countGenome = 16357196
    
        # samtools idxstats LLA073.unjoin_bwa.sort.bam | awk '{s+=$3+$4} END {print s}'
        countStands = 27427625
    
        x = countGenome / (100 * (1-ratio))
    
        # Number of reads required for the standards
        newStands = x * (100 * ratio)
    
        print ('Total number of reads: '  + str(countGenome + newStands))
        print ('Dilution for genome: '    + str(countGenome / (countGenome + newStands)))
        print ('Dilution for standards: ' + str(newStands / (countGenome + newStands)))
        print ('Number of reads required for standards: ' + str(newStands))
    
        genome = 'NA12878_merged.Hg19.sorted.REF_chr21.bam'
        stands = 'LLA073.unjoin_bwa.sort.bam'

        # This is the ratio we'll need to subsample
        ratio = newStands / countStands 

        print ('Subsample ratio: ' + str(ratio))

        # Perform a subsampling (samtools view -s 0.5 -b file.bam > random_half_of_file.bam)
        print ('samtools view -s ' + str(ratio) + ' -b ' + stands + ' > sampled.bam')

        # Run bedtools to perform per-base coverage
        print ('bedtools genomecov -d -ibam ' + genome + ' > genome.txt')
        print ('bedtools genomecov -d -ibam ' + stands + ' > stands.txt')    
        print ('\n')

