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

            file = l1.replace(">","")
            file = file.replace("?","")
            file = file.replace("\n","")
    
            # Create a directory for each sequin
            os.system('mkdir -p ' + (seq_path + file))
        
            w = open(seq_path + file + '/' + file + '.fa', 'w')
            w.write(l1)
            w.write(l2)

def readMixture(file, mix):
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

            tokens = l.strip().split(split)

            #
            # ID  Length MixA  MixB
            #
            
            if (tokens[0] == 'id'):
                continue

            # Create data-structure for the sequin            
            r[tokens[0]] = { 'id':  tokens[0],
                             'len': float(tokens[1]),
                             'A':   float(tokens[2]),
                             'B':   float(tokens[3]),
                           }
    return r

# Generate simulated reads for each sequin for a given mixture
#def simulate(file, basePath, mix='A', min_=0, max_=sys.maxint, c=1, s=60000): # For length normalized 
def simulate(file, basePath, mix, min_=0, max_=sys.maxint, c=1, s=100):
    mixFile = readMixture(file, mix)

    for f in os.listdir(basePath):
        key = f.split('.')[0]

        if key in mixFile:

            # The reads depend on the concentration level and sequin length
            #reads = mixFile[key][mix] / mixFile[key]['len'] # Length normalization
            reads = mixFile[key][mix]
            
            # The concentration needed to be added for the simulation. The number of reads is proportional to the expected concentration.
            reads = c + (s * reads)

            # Length of the sequin
            length = mixFile[key]['len']

            print '\n------------------ ' + key + ' ------------------'
            
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
                print 'Generating: ' + str(reads)

                # Simulate reads from a given sequin
                i  = path + '/' + key + '.fa'
                o1 = path + '/' + key + '.R1.fq'
                o2 = path + '/' + key + '.R2.fq'

                cmd = 'wgsim -s 0 -d 0 -1 100 -2 100 -S ' + str(randint(1,100)) + ' -N ' + str(reads) + ' ' + i + ' ' + o1 + ' ' + o2
                #cmd = 'wgsim -e 0 -r 0 -s 0 -d 0 -1 100 -S ' + str(randint(1,100)) + ' -N ' + str(reads) + ' ' + i + ' ' + o1 + ' /dev/null'

                print cmd
                os.system(cmd)
                    
                # We'll need this command to merge the simulations...
                cmd = 'cp ' + path + '/*.fq ' + basePath
                os.system(cmd)                    
            else:
                print 'Warning: ' + key + ' not generated!'                
        else:
            print '-------- Warning --------: ' + key + ' not found in the mixture!'            

    print('Merging the individual simulations...')
    os.system('cat ' + basePath + '*R1.fq > ' + basePath + 'simulated_1.fastq')
    os.system('cat ' + basePath + '*R2.fq > ' + basePath + 'simulated_2.fastq')

def print_usage():
    print 'Usage: python simulate.py RNA'

if __name__ == '__main__':
    if (len(sys.argv) < 2 or len(sys.argv) > 4):
        print_usage()
    elif (sys.argv[1] == 'RNA'):
        a = ['A1', 'A2', 'A3']              
        b = ['B1', 'B2', 'B3']

        # Simulate replicates for mixture A
        for i in range(0,len(a)):
            split('ATR003.v032.fa', 'RNA_Simulation/')
            simulate('../data/trans/MTR004.v013.csv', 'RNA_Simulation/', 'A')
            os.system('mv RNA_Simulation ' + a[i])

        # Simulate replicates for mixture B
        for i in range(0,len(b)):
            split('ATR003.v032.fa', 'RNA_Simulation/')        
            simulate('../data/trans/MTR004.v013.csv', 'RNA_Simulation/', 'B')
            os.system('mv RNA_Simulation ' + b[i])
    else:
        print_usage()
