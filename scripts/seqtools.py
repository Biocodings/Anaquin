#!/usr/bin/python

#
# This script provides functionality for k-mer
#
#   index -i /Users/tedwong/Desktop/AFU007.v013.index /Users/tedwong/Desktop/AFU007.v013.fa
#   quant -i /Users/tedwong/Desktop/AFU007.v013.index -o /Users/tedwong/Desktop /Users/tedwong/Desktop/A.fq /Users/tedwong/Desktop/B.fq
#

import os
import sys
import vcf
from Bio import SeqIO

class FASTA:
    
    # Name of the chromosome
    name = None
    
    # Base sequenece
    seq = None

class VCF:
    
    # Eg: 7955610
    pos = None
    
    # Reference allele
    ref = None
    
    # Alternative alleles
    alt = None
    
    # Is this a SNP?
    isSNP = None
    
    # Is this an indel?
    isIndel = None

def readFA(file):
    parser = SeqIO.parse(open(file), 'fasta')
    x = {}

    for fasta in parser:
        name, seq = fasta.id, str(fasta.seq)
        
        f = FASTA()
        
        f.name = name
        f.seq  = seq

        x[f.name] = f
        
    return x

def readVCF(file):
    reader = vcf.Reader(open(file, 'r'))
    r = {}

    for record in reader:
        
        x = VCF()
        
        x.pos      = record.POS
        x.ref      = record.REF
        x.alt      = record.ALT
        x.isSNP    = record.is_snp
        x.is_indel = record.is_indel
        
        r[str(x.pos)] = x

    return r

#
# Usage: python seqtools reverse <genome.fa> <variants.vcf>
#
if __name__ == '__main__':

    # Eg: genome.fa
    fFile = '/Users/tedwong/Desktop/chrF.fa'  #sys.argv[2]

    # Eg: variants.vcf
    vFile = '/Users/tedwong/Sources/QA/data/VarQuin/AVA009.v032.vcf' #sys.argv[3]

    fData = readFA(fFile)
    vData = readVCF(vFile)

    











    
