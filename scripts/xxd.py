#
# This script converts resources into C++ headers such that the software can compile with.
#

#!/usr/bin/python

import os
import sys
import subprocess

os.system('rm -rf src/resources/*')

def xxd(src, dst):
    os.system('xxd -i ' + src + ' ' + ' ' + dst)

data  = [ 'data/manual.txt',    
          'data/linear.R',

           # ---------- META ----------

          'data/meta/m_manual.txt',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v6.mix.csv',
          'data/meta/META.v1.tab.bed',

           # ---------- RNA ----------

          'data/rna/r_manual.txt',
          'data/rna/RNA.v0.chrT.fa',  # Reference for the chromosome
          'data/rna/RNA.v1.chrT.bed', # Reference BED for the chromosome
          'data/rna/RNA.v1.chrT.gtf', # Reference GTF for the chromosome
          'data/rna/RNA.v1.mix.csv',  # Mixture CSV file
          'data/rna/RNA.v1.tab.fa',   # Reference for the sequins

           # ---------- DNA ----------

          'data/dna/d_manual.txt',
          'data/dna/DNA.v0.chrT.fa', # Reference for the chromosome
          'data/dna/DNA.v0.vcf',     # Reference structural variations
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
