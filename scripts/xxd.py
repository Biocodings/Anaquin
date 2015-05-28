#
# This script converts resources into C++ headers such that the software can compile with.
#

#!/usr/bin/python

import os
import sys
import subprocess

os.system('rm -rf src/data/*')

def xxd(src, dst):
    os.system('xxd -i ' + src + ' ' + ' ' + dst)

docs  = [ 'docs/manual.txt' ]
data  = [ 'data/silico.fa',

           # ---------- META ----------

          'data/meta/META.tab.fa',
          'data/meta/META.mix.csv',

           # ---------- RNA ----------

          'data/rna/RNA.tab.fa',
          'data/rna/RNA.ref.gtf',
          'data/rna/RNA.ref.bed',
          'data/rna/RNA.mix.csv',

           # ---------- DNA ----------

          'data/dna/DNA.mix.csv',
          'data/dna/DNA.var.vcf',          
          'data/dna/DNA.ref.bed',
          'data/dna/DNA.tab.fa',
        ]
tests = [ ]

r = data + docs
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/data/' + file)
