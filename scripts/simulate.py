#!/usr/bin/python

import os
import sys
import subprocess

def d1_path():
    return '../data/d1/'

def d1_seq_path():
    return d1_path() + 'seqs/'

def d1_read_path():
    return d1_path() + 'reads/'

def d1_rna_mix_a_path():
    return '../data/RNA/Standard_A.csv'

def d_sequins():
    return '../data/DNA/DNA.tab.fa'

def r_sequins():
    return '../data/RNA/RNA.tab.fa'

def d_standards():
    return '../data/DNA/DNA_Standards_Analysis.txt'

def r_standards():
    return '../data/RNA/RNA_Standards_Analysis.txt'

# Split a file of sequin into individual sequins
def split_sequins(file):
    os.system('mkdir -p ' + d1_seq_path())

    with open(file) as f:
        while True:
            l1 = f.readline()
            l2 = f.readline()
            
            if (not l2):
                break
            
            file = l1.replace(">","")
            file = file.replace("?","")
            file = file.replace("\n","")
            
            w = open(d1_seq_path() + file + '.fa', 'w')
            w.write(l1)
            w.write(l2)

def csv_to_stand(ts):
    x = {
        'gene': ts[2],
        'ts': ts[3],
        'amount': float(ts[5]),
        'ratio': ts[4]
    }

    return x

def read_standards(file):
    ps = {}
    with open(file) as f:
        l = f.readline()
        while True:
            l = f.readline()
            if (not l):
                break
            
            tokens = l.split('\t')
            ps[tokens[1]] = { 'id': tokens[1],
                              'con_a': float(tokens[3]),
                              'con_b': float(tokens[4]) }

    return ps

# Generate simulated reads for each sequin
def simulate_reads(file):
    os.system('mkdir -p ' + d1_read_path())
    ps = read_standards(file)

    for f in os.listdir(d1_seq_path()):
        # Name of the transcript
        ts = f.split('.')[0]
        
        if ts in ps:
            amount = ps[ts]['con_a']
            
            # Just to make sure ...
            n = amount + 10000
            
            print '------------------ ' + ts + ' ------------------'
            print 'Generating: ' + str(n)

            # Command: [wgsim -d 400 -N 10000 -1 101 -2 101 ${X} ${X}.R1.fq ${X}.R2.fq]
            
            i  = d1_seq_path()  + ts + '.fa'
            o1 = d1_read_path() + ts + '.R1.fq'
            o2 = d1_read_path() + ts + '.R2.fq'

            # Simulate reads from a given sequin
            cmd = 'wgsim -d 400 -N ' + str(n) + ' -1 101 -2 101 ' + i + ' ' + o1 + ' ' + o2

            os.system(cmd)

    print('Merging the individual simulations...')
    os.system('cat ' + d1_seq_path() + '*.fa > ' + d1_seq_path() + 'simulated.fq')

if __name__ == '__main__':
    if (len(sys.argv) != 2):
        print 'Usage: python simulate.py RNA|DNA'
    elif (sys.argv[1] == 'DNA'):
        split_sequins(d_sequins())
        simulate_reads(d_standards())
    elif (sys.argv[1] == 'RNA'):
        split_sequins(r_sequins())
        simulate_reads(r_standards())
    else:
        print 'Usage: python simulate.py RNA|DNA'





