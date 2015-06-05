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

docs  = [ 'resources/manual.txt' ]
data  = [ 'resources/silico.fa',

           # ---------- META ----------

          'resources/meta/META.tab.fa',
          'resources/meta/META.mix.csv',
          'resources/meta/META.ref.bed',          

           # ---------- RNA ----------

          'resources/rna/RNA.tab.fa',
          'resources/rna/RNA.ref.gtf',
          'resources/rna/RNA.ref.bed',
          'resources/rna/RNA.mix.csv',

           # ---------- DNA ----------

          'resources/dna/DNA.mix.csv',
          'resources/dna/DNA.var.vcf',          
          'resources/dna/DNA.ref.bed',
          'resources/dna/DNA.tab.fa',
        ]
tests = [ ]

r = data + docs
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/data/' + file)
