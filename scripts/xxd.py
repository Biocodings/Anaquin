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
          'data/chrT.v1.fa',

           # ---------- RNA ----------

          'data/rna/RNA.usage.txt',
          'data/rna/RNA.v1.bed',
          'data/rna/RNA.v1.gtf',
          'data/rna/RNA.v4.1.mix',
          'data/rna/RNA.v1.2.tab.fa',

           # ---------- META ----------

          'data/meta/META.usage.txt',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v6.mix.csv',
          'data/meta/META.v1.tab.bed',

           # ---------- DNA ----------

          'data/dna/DNA.usage.txt',
          'data/dna/DNA.v3.mix.csv',
          'data/dna/DNA.variant.bed',
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
