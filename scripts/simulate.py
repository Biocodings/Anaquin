#!/usr/bin/python

import os
import subprocess

def d1_path():
    return '../data/d1/'

def d1_seq_path():
    return d1_path() + 'seqs/'

def d1_read_path():
    return d1_path() + 'reads/'

def d1_rna_mix_a_path():
    return '../data/RNA/Standard_A.csv'

# Split a file of sequin into individual sequins
def split_sequins():
    os.system('mkdir -p ' + d1_seq_path())

    with open('../data/RNA/RNAsequins.fa') as f:
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
        while True:
            l1 = f.readline()
            l2 = f.readline()
            
            if (not l1 or not l2):
                break
            
            t1 = csv_to_stand(l1.split(','))
            t2 = csv_to_stand(l2.split(','))
            
            ps[t1['ts']] = t1
            ps[t2['ts']] = t2

    return ps

# Generate simulated reads for each sequin
def simulate_reads(file):
    os.system('mkdir -p ' + d1_read_path())

    ps = read_standards(file)
    
    for f in os.listdir(d1_seq_path()):
        # Name of the transcript
        ts = f.split('.')[0]
        
        print ts + ' found'

        if ts in ps:
            amount = ps[ts]['amount']
            
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

if __name__ == '__main__':
    split_sequins()
    simulate_reads(d1_rna_mix_a_path())
