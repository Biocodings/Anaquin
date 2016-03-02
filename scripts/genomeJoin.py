#!/usr/bin/python

import os
import random

if __name__ == '__main__':

    names = []
    seqs  = []
    
    with open("CLA004.v013.fa", "r") as f:
        for line in f:
            if (line[0] == '>'):
                line = line.replace('>', '')                
                names.append(line.strip())
            else:
                seqs.append(line.strip())

    writer = open('joined.fa', 'w')
    
    writer.write('>chrT\n')
    
    for seq in seqs:    
        writer.write(seq)
        
    # Run: fasta_formatter -i joined.fa -o Temp.fa -w 60; mv Temp.fa joined.fa
        
    writer.write('\n')
    writer.close()