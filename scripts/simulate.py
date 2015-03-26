#!/usr/bin/python

import os
import sys
import math
import subprocess

def dna_path():
    return '../data/dna_sim/'

def dna_seq_path():
    return dna_path() + 'seqs/'

def dna_read_path():
    return dna_path() + 'reads/'

def rna_path():
    return '../data/rna_sim/'

def rna_seq_path():
    return rna_path() + 'seqs/'

def rna_read_path():
    return rna_path() + 'reads/'

def r_sequins():
    return '../data/RNA/RNA.tab.fa'

def d_sequins():
    return '../data/DNA/DNA.tab.fa'

def d_standards():
    return '../data/DNA/DNA_Standards_Analysis.txt'

def r_standards():
    return '../data/RNA/RNA_Standards_Analysis.txt'

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

def read_standards(file):
    ps = {}
    with open(file) as f:
        l = f.readline()
        while True:
            l = f.readline()
            if (not l):
                break
            
            tokens = l.split('\t')
            ps[tokens[1]] = { 'id':    tokens[1],
                              'group': tokens[2],
                              'con_a': float(tokens[3]),
                              'con_b': float(tokens[4]),
                              'ratio': float(tokens[5]),
                              'logr' : float(tokens[6]) }

    return ps

# Generate simulated reads for each sequin
def simulate_reads(file, seq_path, read_path):
    os.system('mkdir -p ' + read_path)
    ps = read_standards(file)

    for f in os.listdir(seq_path):
        ts = f.split('.')[0]
        
        if ts in ps:
            na = ps[ts]['con_a']
            nb = ps[ts]['con_b']
            ratio = math.log(na / nb, 2)
            
            if (math.fabs(ratio - ps[ts]['logr']) > 0.5):
                raise Exception('Inconsistence mixture ratio: ' + ps[ts]['id'])

            # Multiply the concentration by a constant
            na = (1000 * na) + 1000
            
            print '------------------ ' + ts + ' ------------------'
            print 'Generating: ' + str(na)
            
            # Command: wgsim -d 400 -N 10000 -1 101 -2 101 ${X} ${X}.R1.fq ${X}.R2.fq
            
            i  = seq_path  + ts + '.fa'
            o1 = read_path + ts + '.R1.fq'
            o2 = read_path + ts + '.R2.fq'

            # Simulate reads from a given sequin
            cmd = 'wgsim -d 400 -N ' + str(int(na)) + ' -1 101 -2 101 ' + i + ' ' + o1 + ' ' + o2

            print cmd
            os.system(cmd)

    print('Merging the individual simulations...')
    os.system('cat ' + seq_path + '*.fa > ' + read_path + 'simulated.fq')

if __name__ == '__main__':
    if (len(sys.argv) != 2):
        print 'Usage: python simulate.py RNA|DNA'
    elif (sys.argv[1] == 'DNA'):
        split_sequins(d_sequins(), dna_seq_path())
        simulate_reads(d_standards(), dna_seq_path(), dna_read_path())
    elif (sys.argv[1] == 'RNA'):
        split_sequins(r_sequins(), rna_seq_path())
        simulate_reads(r_standards(), rna_seq_path(), rna_read_path())
    else:
        print 'Usage: python simulate.py RNA|DNA'

    # Reads have been simulated and written to simulated.fq.
    #
    # Build an index for bowtie:
    #    --> bowtie2-build -f ../chromo/ChrT.5.10.fa index
    #
    # Example workflow for RNA:
    #
    #    1. tophat2 -p 8 -G ../chromo/chromo.gtf -o aligned index seqs/simulated.fq
    #
    # Example workflow for DNA:
    #
    #    1. tophat2 -p 8 -G ../chromo/chromo.gtf -o aligned index seqs/simulated.fq
    #