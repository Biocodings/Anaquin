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

docs  = [ 'docs/manual.txt' ]
data  = [ 'data/silico.fa',

          'scripts/linear.R',

           # ---------- META ----------

          'data/meta/META.v1.ref.fa',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v6.mix.csv',
          'data/meta/META.v1.ref.bed',

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
    xxd(r[i], 'src/resources/' + file)
