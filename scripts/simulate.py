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
            # ID  Length Mix_A  Mix_B
            #
            
            if (tokens[0] == 'id'):
                continue
                
            # Create data-structure for the sequin            
            r[tokens[0]] = { 'id': tokens[0],                            
                             'a':  float(tokens[2]),
                             'b':  float(tokens[3]),
                           }
    return r

# Generate simulated reads for each sequin for a given mixture
def simulate(file, basePath, mix='A', c=0, s=1, tool='wgsim'):
    mix = readMixture(file, mix)

    for f in os.listdir(basePath):
        key = f.split('.')[0]

        if key in mix:            
            # The concentration level depends on the level
            if mix == 'A':
                con = mix[key]['a']
            else:
                con = mix[key]['b']

            # The concentration needed to be added for the simulation
            con = c + (s * con)

            print '\n------------------ ' + key + ' ------------------'
            
            # Path for the sequin
            path  = basePath  + key

            # This is the number of reads that we'll need
            con = int(con)

            # Don't bother if the abundance is too low or too high
            if (con > 1 and con < 200000):
                print 'Generating: ' + str(con)

                # Simulate reads from a given sequin
                if tool == 'wgsim':
                    i  = path + '/' + key + '.fa'
                    o1 = path + '/' + key + '.R1.fastq'
                    o2 = path + '/' + key + '.R2.fastq'                    

                    # Command: wgsim -e 0 -N 5151 -1 101 -2 101 ${X} ${X}.R1.fastq ${X}.R2.fastq            
                    cmd = 'wgsim -r 0 -s 0 -S ' + str(randint(1,100)) + ' -N ' + str(con) + ' ' + i + ' ' + o1 + ' ' + o2

                    print cmd
                    os.system(cmd)
                    
                    # We'll need this command to merge the simulations...
                    cmd = 'cp ' + path + '/*.fastq ' + basePath
                    os.system(cmd)                    
                else:                    
                    os.system('mkdir -p ' + os.getcwd() + '/Sherman/' + key)
                    #cmd = '~/scripts/Sherman -cr 0 -e 0 -n ' + str(con) + ' -l 101 --genome_folder ' + os.getcwd() + '/Sherman/' + key
                    cmd = '../Sherman -cr 0 -e 0 -n ' + str(con) + ' -l 101 --genome_folder ' + os.getcwd() + '/Sherman/' + key

                    print cmd
                    os.system(cmd)

                    # Move the generated FASTQ to the proper directory
                    print ('mv simulated.fastq ' + os.getcwd() + '/Sherman/' + key)
                    os.system('mv simulated.fastq ' + os.getcwd() + '/Sherman/' + key + "/")
                    
                    os.system('find Sherman/ -name *.fastq -exec cat {} + > simulated.fastq')
            else:
                print 'Warning: ' + key + ' not generated!'                
        else:
            print '-------- Warning --------: ' + key + ' not found in the mixture!'            

    if (tool == 'wgsim'):
        print('Merging the individual simulations...')
        os.system('cat ' + basePath + '*R1.fastq > ' + basePath + 'simulated_1.fastq')
        os.system('cat ' + basePath + '*R2.fastq > ' + basePath + 'simulated_2.fastq')
        #os.system('rm '  + basePath + '/*R1.fastq')
        #os.system('rm '  + basePath + '/*R2.fastq')

def print_usage():
    print 'Usage: python simulate.py RNA|DNA|META|RNAFlat'

if __name__ == '__main__':
    if (len(sys.argv) < 2 or len(sys.argv) > 4):
        print_usage()
    elif (sys.argv[1] == 'RNAFlat'):
        a = ['RNA_A_1', 'RNA_A_2', 'RNA_A_3']              
        b = ['RNA_B_1', 'RNA_B_2', 'RNA_B_3']

        for i in range(0,len(a)):
            split('../data/rna/RNA.v1.fa', 'RNA_Simulation/')
            simulate('../tests/data/rna/flat.mix.csv', 'RNA_Simulation/', 'A')
            os.system('mv RNA_Simulation ' + a[i])

        for i in range(0,len(b)):
            split('../data/rna/RNA.v1.fa', 'RNA_Simulation/')        
            simulate('../tests/data/rna/flat.mix.csv', 'RNA_Simulation/', 'B')
            os.system('mv RNA_Simulation ' + b[i])
    elif (sys.argv[1] == 'RNA'):
        a = ['RNA_A_1', 'RNA_A_2', 'RNA_A_3']              
        b = ['RNA_B_1', 'RNA_B_2', 'RNA_B_3']

        for i in range(0,len(a)):
            split('../data/rna/RNA.v1.fa', 'RNA_Simulation/')
            simulate('../data/rna/RNA.v4.1.mix', 'RNA_Simulation/', 'A')
            os.system('mv RNA_Simulation ' + a[i])

        for i in range(0,len(b)):
            split('../data/rna/RNA.v1.fa', 'RNA_Simulation/')        
            simulate('../data/rna/RNA.v4.1.mix', 'RNA_Simulation/', 'B')
            os.system('mv RNA_Simulation ' + b[i])
    elif (sys.argv[1] == 'META'):
        split('../data/meta/META.v1.tab.fa', 'META_A/')
        split('../data/meta/META.v1.tab.fa', 'META_B/')

        # Generate simulation for mixture A (5% of the origianl concentration to save time)
        simulate('../data/meta/META.v6.mix.csv', 'META_A/', 'A', 0, 1.0)

        # Generate simulation for mixture B (5% of the origianl concentration to save time)
        simulate('../data/meta/META.v6.mix.csv', 'META_B/', 'B', 0, 1.0)
    else:
        print_usage()