#!/usr/bin/python

import os
import sys
import math
import subprocess

def dna_path():
    return 'DNA_Simulation/'

def rna_path():
    return 'RNA_Simulation/'

def seq_path():
    return rna_path() + 'seqs/'

def read_path():
    return rna_path() + 'reads/'

def r_sequins():
    return '../data/rna/RNA.tab.fa'

def d_sequins():
    return '../data/dna/DNA.tab.fa'

def r_mixtures():
    return '../data/rna/rna_mixtures.csv'

# Split a file of sequin into individual sequins
def split_sequins(file, seq_path):
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

def read_mixture(file, mix):
    r = {}
    with open(file) as f:
        l = f.readline()
        while True:
            l = f.readline()
            if (not l):
                break

            tokens = l.strip().split(',')
            
            #
            # Ref, Var, Group, RLen, VLen, AimA, AimB, RatioR, RatioV
            #

            ratio_r = float(tokens[7]);
            ratio_v = float(tokens[8]);
            ratio_t = ratio_r + ratio_v;

            aim   = float(tokens[5] if mix == 'A' else tokens[6])
            ratio = float(ratio_r   if mix == 'A' else ratio_v)      
            r_raw = float(aim * (ratio / ratio_t))
            v_raw = float(aim - r_raw)

            r_fpkm = r_raw / (1000 / float(tokens[3]))
            v_fpkm = -1 if tokens[4] == '' else v_raw / (1000 / float(tokens[4]))

            # Create data-structure for the reference sequin
            r[tokens[0]] = { 'id':    tokens[0],
                             'group': tokens[2],
                             'raw':   r_raw,
                             'fpkm':  r_fpkm,
                           }

            # Create data-structure for the variant sequin
            r[tokens[1]] = { 'id':    tokens[1],
                             'group': tokens[2],
                             'raw':   v_raw,
                             'fpkm':  v_fpkm,
                           }
    return r

# Generate simulated reads for each sequin for a given mixture
def simulate_reads(file, seq_path, read_path, min_, max_, mix):
    os.system('mkdir -p ' + read_path)
    d = read_mixture(file, mix)

    for f in os.listdir(seq_path):
        ts = f.split('.')[0]
    
        if ts in d:
            fpkm = d[ts]['fpkm']
            
            # This can be changed to simulate sequencing depth

            s = 1
            c = 0

            # Multiply the concentration by a constant (applies to all sequins)
            fpkm = c + (s * fpkm)

            na = max(min_, fpkm)
            na = min(max_, fpkm)
            
            print '\n------------------ ' + ts + ' ------------------'
            
            # Command: wgsim -d 400 -N 5151 -1 101 -2 101 ${X} ${X}.R1.fastq ${X}.R2.fastq
            
            i  = seq_path  + ts + '.fa'
            o1 = read_path + ts + '.R1.fastq'
            o2 = read_path + ts + '.R2.fastq'

            # Don't bother if the abundance is too low
            if (int(na) > 5):
                print 'Generating: ' + str(fpkm)
                
                # Simulate reads from a given sequin
                cmd = 'wgsims -d 400 -N ' + str(int(fpkm)) + ' -1 101 -2 101 ' + i + ' ' + o1 + ' ' + o2

                print cmd
                os.system(cmd)
        else:
            print 'Ignore ' + ts + ' (not found)

    print('Merging the individual simulations...')
    os.system('cat ' + read_path + '*R1.fastq > ' + read_path + 'simulated_1.fastq')
    os.system('cat ' + read_path + '*R2.fastq > ' + read_path + 'simulated_2.fastq')

def print_usage():
    print 'Usage: python simulate.py RNA_A|RNA_B'

if __name__ == '__main__':
    if (len(sys.argv) != 2):
        print_usage()
    elif (sys.argv[1] == 'RNA_A'):
        split_sequins(r_sequins(), seq_path())
        simulate_reads(r_mixtures(), seq_path(), read_path(), 0, sys.maxint, 'A')
    elif (sys.argv[1] == 'RNA_B'):
        split_sequins(r_sequins(), seq_path())
        simulate_reads(r_mixtures(), seq_path(), read_path(), 0, sys.maxint, 'B')
    #elif (sys.argv[1] == 'DNA_T'):
    #    split_sequins(d_sequins(), dna_seq_path())
    #    simulate_reads(d_standards(), dna_seq_path(), dna_read_path(), 120000, 120000)
    #elif (sys.argv[1] == 'DNA'):
    #    split_sequins(d_sequins(), dna_seq_path())
    #    simulate_reads(d_standards(), dna_seq_path(), dna_read_path(), 0, sys.maxint)
    else:
        print_usage()

    # Reads have been simulated and written to simulated.fq.
    #
    # Build an index for bowtie:
    #    --> bowtie2-build -f ../chromo/ChrT.5.10.fa index
    #
    # Example workflow for RNA:
    #
    #    1. tophat2 -p 8 -G ../chromo/chromo.gtf -o aligned index seqs/simulated.fq
