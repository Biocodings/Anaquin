#!/usr/bin/python

#
# This script converts VarScan output file (SNP) to BED. The script is designed solely for quick sequins
# validation.
#

import os
import sys

#
# BED:
#
#   <chrT>  <Position>  <Position>  <ID>
#

#
# VarScan:
#
#   <Chrom> <Position>  <Ref>   <Cons>
#

starts = []
names  = []

with open('tests/data/varscan.tab') as f:
    while True:
        l = f.readline()
        if (not l):
            break
            
        toks = l.split()
        
        if (toks[0] == 'chrT'):
            if (len(toks) >= 19):
                if (toks[2] != toks[18]):
                    starts.append(toks[1])
                    names.append(toks[2] + '_' + toks[18])

out = open('Ira.bed', 'w')      
for i in range(len(starts)):
    out.write('chrT\t' + str(starts[i]) + '\t' + str(starts[i]) + '\t' + names[i] + '\n')

out.close()

