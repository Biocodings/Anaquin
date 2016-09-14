#!/usr/bin/python

import os
import sys
import math
import subprocess
from random import randint

# Split a file of sequin into individual sequins
def split(file, seq_path):
    os.system('mkdir -p ' + seq_path)

    with open(file) as f:
        while True:
            l1 = f.readline()
            l2 = f.readline()            
            if (not l2):
                break

            file = l1.replace(">", "")
            file = file.replace("?", "")
            file = file.replace("\n", "")
    
            # Create a directory for each sequin
            os.system('mkdir -p ' + (seq_path + file))
        
            w = open(seq_path + file + '/' + file + '.fa', 'w')
            w.write(l1)
            w.write(l2)

def readMixture(file):
    r = {}

    with open(file) as f:
        l = f.readline()

        if (',' in l):
            split = ','
        else:
            split = '\t'

        while True:
            l = f.readline()
            if (not l):
                break

            toks = l.strip().split(split)

            #
            # ID  Length MixA  MixB
            #

            if (toks[0] == 'id'):
                continue

            # Create data-structure for the sequin            
            r[toks[0]] = { 'id':  toks[0],
                           'len': float(toks[1]),
                           'A':   float(toks[2]),
                           'B':   float(toks[3]),
                         }
    return r

# Generate simulated reads for each sequin for a given mixture
def simulate(file, basePath, mix, min_=0, max_=sys.maxint, c=0, s=1000, isTooLowError=True):
    mixFile = readMixture(file)
    covs = []

    for f in os.listdir(basePath):
        key = f.split('.')[0]

        if key in mixFile:
            # Length of the sequin
            length = mixFile[key]['len']

            # The reads depend on the concentration level and sequin length
            #reads = mixFile[key][mix] / mixFile[key]['len'] # Length normalization
            reads = mixFile[key][mix]

            # The concentration needed to be added for the simulation. The number of reads is proportional to the expected concentration.
            reads = c + (s * reads)

            #print '\n------------------ ' + key + ' ------------------'
            
            # Path for the sequin
            path  = basePath  + key

            # This is the number of reads that we'll need
            reads = int(reads)

            #
            # The number of reads need to be adjusted for the sequin length.
            # The implementation borrows from Wendy's script. For example,
            #
            #   - length is 1689
            #   - reads is 468750
            #
            # We would calculate: 468750 * (1689 / 1000) to get reads per KB.
            #
            
            reads = max(min_, reads)
            reads = min(max_, reads)
            reads = reads * (length / 1000)

            # Don't bother if the abundance is too low or too high
            if (reads > 1):
                #print 'Generating: ' + str(reads)

                # Simulate reads from a given sequin
                i  = path + '/' + key + '.fa'
                o1 = path + '/' + key + '.R1.fq'
                o2 = path + '/' + key + '.R2.fq'

                #cmd = 'wgsim -1150 -2150 -S ' + str(randint(1,100)) + ' -N ' + str(reads) + ' ' + i + ' ' + o1 + ' ' + o2 + ' > /dev/null'
                cmd = 'wgsim -1150 -2150 -d 130 -S ' + str(randint(1,100)) + ' -N ' + str(reads) + ' ' + i + ' ' + o1 + ' ' + o2 + ' > /dev/null'
                                
                # What's the estimated coverage?
                cov = (150 * 2 * reads) / length
                
                covs.append(cov)
                
                #print cmd
                os.system(cmd)
                
                if (os.stat(o1).st_size == 0):
                    raise Exception(key + ' not generated!')
                if (os.stat(o2).st_size == 0):
                    raise Exception(key + ' not generated!')
                    
                # We'll need this command to merge the simulations...
                cmd = 'cp ' + path + '/*.fq ' + basePath
                os.system(cmd)                    
            else:
                if isTooLowError:
                    raise Exception('Error: ' + key + ' not generated! Reads: ' + str(reads))
                else:
                    print 'Warning: ' + key + ' not generated!'
        else:
            print '-------- Warning --------: ' + key + ' not found in the mixture!'            

    #print('Merging the individual simulations...')
    os.system('cat ' + basePath + '*R1.fq > ' + basePath + 'simulated_1.fastq')
    os.system('cat ' + basePath + '*R2.fq > ' + basePath + 'simulated_2.fastq')

    print('Minimum coverage for mixture ' + str(mix) + ': ' + str(min(covs)))
    print('Maximum coverage for mixture ' + str(mix) + ': ' + str(max(covs)))

def print_usage():
    print 'Usage: python simulate.py RNA'

#
# Eg: python simulate.py RNA A.R.2.fa NewMixAB.csv
#
#   data <- read.csv('/Users/tedwong/Desktop/ABCD.csv', sep=',')
#   write.table(data, sep='\t', quote=FALSE, file='/Users/tedwong/Desktop/MixAB.csv', row.names=FALSE)
#
#   bowtie2-build A.R.4.fa A.R.4
#
#   tophat2 -o A1 -p 6 A.R.4 A1/simulated_1.fastq A1/simulated_2.fastq
#   tophat2 -o A2 -p 6 A.R.4 A2/simulated_1.fastq A2/simulated_2.fastq
#   tophat2 -o A3 -p 6 A.R.4 A3/simulated_1.fastq A3/simulated_2.fastq
#   tophat2 -o B1 -p 6 A.R.4 B1/simulated_1.fastq B1/simulated_2.fastq
#   tophat2 -o B2 -p 6 A.R.4 B2/simulated_1.fastq B2/simulated_2.fastq
#   tophat2 -o B3 -p 6 A.R.4 B3/simulated_1.fastq B3/simulated_2.fastq
#
#   cufflinks -p 15 -G ARN020_v001.gtf -o A1/G A1/accepted_hits.bam
#   cufflinks -p 15 -G ARN020_v001.gtf -o A2/G A2/accepted_hits.bam
#   cufflinks -p 15 -G ARN020_v001.gtf -o A3/G A3/accepted_hits.bam
#   cufflinks -p 15 -G ARN020_v001.gtf -o B1/G B1/accepted_hits.bam
#   cufflinks -p 15 -G ARN020_v001.gtf -o B2/G B2/accepted_hits.bam
#   cufflinks -p 15 -G ARN020_v001.gtf -o B3/G B3/accepted_hits.bam
#
if __name__ == '__main__':    
    if (len(sys.argv) < 2 or len(sys.argv) > 4):
        print_usage()
    elif (sys.argv[1] == 'RNA'):
        fasta = sys.argv[2]
        mixAB = sys.argv[3]
    
        print('FASTA: ' + fasta)
        print('Mixture: ' + mixAB)
        
        # We'll simulate three replicates for mixture A
        a = ['A1', 'A2', 'A3']
        
        # We'll also simulate three replicates for mixture B
        b = ['B1', 'B2', 'B3']
        
        os.system('rm -rf A1 A2 A3 B1 B2 B3 RNA_Simulation')

        # Simulate replicates for mixture A
        for i in range(0,len(a)):
            split(fasta, 'RNA_Simulation/')
            simulate(mixAB, 'RNA_Simulation/', 'A')
            os.system('mv RNA_Simulation ' + a[i])

        # Simulate replicates for mixture B
        for i in range(0,len(b)):
            split(fasta, 'RNA_Simulation/')        
            simulate(mixAB, 'RNA_Simulation/', 'B')
            os.system('mv RNA_Simulation ' + b[i])
    else:
        print_usage()
