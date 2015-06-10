#!/usr/bin/python

import os
import sys
import math
import subprocess
from random import randint

def rna_path():
    return 'RNA_Simulation/'

def seq_path(base):
    return base + 'seqs/'

def read_path(base):
    return base + 'reads/'

def r_sequins():
    return '../data/rna/RNA.tab.fa'

def d_sequins():
    return '../data/dna/DNA.tab.fa'

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
    
            w = open(seq_path + file + '.fa', 'w')
            w.write(l1)
            w.write(l2)

def readMixture(file, mix):
    r = {}

    with open(file) as f:
        l = f.readline()

        while True:
            l = f.readline()
            if (not l):
                break

            tokens = l.strip().split(',')

            #
            # ID  Mix_A  Mix_B
            #
            
            if (tokens[0] == 'id'):
                continue

            # Create data-structure for the sequin            
            r[tokens[0]] = { 'id': tokens[0],
                             'a':  float(tokens[1]),
                             'b':  float(tokens[2]),
                           }
    return r

# Generate simulated reads for each sequin for a given mixture
def simulate(file, seq_path, read_path, mix, tool='wgsim', c=0, s=1):
    os.system('mkdir -p ' + read_path)
    mix = readMixture(file, mix)

    for f in os.listdir(seq_path):
        key = f.split('.')[0]

        if key in mix:
            if mix == 'A':
                con = mix[key]['a']
            else:
                con = mix[key]['b']

            # The concentration needed to be added for the simulation
            con = c + (s * con)

            print '\n------------------ ' + key + ' ------------------'
            
            # Command: wgsim -e 0 -d 400 -N 5151 -1 101 -2 101 ${X} ${X}.R1.fastq ${X}.R2.fastq
            
            i  = seq_path  + key + '.fa'
            o1 = read_path + key + '.R1.fastq'
            o2 = read_path + key + '.R2.fastq'

            # This is the number of reads that we'll need
            con = int(con)

            # Don't bother if the abundance is too low
            if (con > 1):
                print 'Generating: ' + str(con)

                # Simulate reads from a given sequin
                if tool == 'wgsim':
                    cmd = 'wgsim -r 0 -S ' + str(randint(1,100)) + '  -d 400 -N ' + str(int(con)) + ' -1 101 -2 101 ' + i + ' ' + o1 + ' ' + o2
                else:
                    os.system('mkdir -p ' + os.getcwd() + '/Sherman/' + key)
                    cmd = '~/scripts/Sherman -cr 0 -e 0 -n ' + str(con) + ' -l 101 --genome_folder ' + os.getcwd() + '/Sherman/' + key

                print cmd
                os.system(cmd)
            else:
                print '-------- Warning --------: ' + key + ' not generated!'                
        else:
            print '-------- Warning --------: ' + key + ' not found!'

    print('Merging the individual simulations...')
    os.system('cat ' + read_path + '*R1.fastq > ' + read_path + 'simulated_1.fastq')
    os.system('cat ' + read_path + '*R2.fastq > ' + read_path + 'simulated_2.fastq')

def print_usage():
    print 'Usage: python simulate.py RNA|DNA|META|Sherman'

if __name__ == '__main__':
    if (len(sys.argv) < 2 or len(sys.argv) > 4):
        print_usage()
    elif (sys.argv[1] == 'Sherman'):
        split(sys.argv[2], 'Sherman/seqs/')
        simulate(sys.argv[3], 'Sherman/seqs/', 'Sherman/reads/', 'A', 0 ,1)
        #simulate(sys.argv[3], 'Sherman/seqs/', 'Sherman/reads/', 'B', 0 ,1)                
    elif (sys.argv[1] == 'RNA'):
        a = ['RNA_A_1_1', 'RNA_A_1_2', 'RNA_A_1_3']              
        for i in range(0,len(a)):
            split_sequins(r_sequins(), seq_path(rna_path()))        
            simulate_reads(r_mixtures(), seq_path(rna_path()), read_path(rna_path()), 0, sys.maxint, 'A')
            os.system('mv RNA_Simulation ' + a[i])
        b = ['RNA_B_10_1', 'RNA_B_10_2', 'RNA_B_10_3']
        for i in range(0,len(b)):
            split_sequins(r_sequins(), seq_path(rna_path()))        
            simulate_reads(r_mixtures(), seq_path(rna_path()), read_path(rna_path()), 0, sys.maxint, 'B', 0, 10)
            os.system('mv RNA_Simulation ' + b[i])
    elif (sys.argv[1] == 'DNA'):
        pass
    elif (sys.argv[1] == 'META'):
        split('../data/meta/META.ref.fa', 'META_Simulation_A/seqs/')
        split('../data/meta/META.ref.fa', 'META_Simulation_B/seqs/')
        
        # Generate simulation for mixture A (5% of the origianl concentration to save time)
        simulate('../data/meta/META.mix.csv', 'META_Simulation_A/seqs/', 'META_Simulation_A/reads/', 'A', 0, 0.05)
        
        # Generate simulation for mixture B (5% of the origianl concentration to save time)
        simulate('../data/meta/META.mix.csv', 'META_Simulation_B/seqs/', 'META_Simulation_B/reads/', 'B', 0, 0.05)

        os.system('velveth A 31 -fastq -shortPaired META_Simulation_A/reads/simulated_1.fastq META_Simulation_A/reads/simulated_2.fastq')
        #os.system('velvetg A')
        os.system('velvetg A -exp_cov auto')
        #os.system('meta-velvetg A')
        os.system('blat ../data/meta/META.ref.fa A/contigs.fa A/align.psl')

        os.system('velveth B 21 -fastq -shortPaired META_Simulation_B/reads/simulated_1.fastq META_Simulation_B/reads/simulated_2.fastq')    
        #os.system('velvetg B')
        os.system('velvetg B -exp_cov auto')
        #os.system('meta-velvetg A')        
        os.system('blat ../data/meta/META.ref.fa B/contigs.fa B/align.psl')
    else:
        print_usage()