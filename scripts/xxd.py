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

           # ------- Fusion -------

           'data/fusion/FUS.v3.csv',
           'data/fusion/FUS.v1.ref',

           # ------- Ladder -------

           'data/ladder/CON.v3.mix.csv',

           # ---------- RNA ----------

          'data/rna/RNA.usage.txt',
          'data/rna/RNA.v1.gtf',
          'data/rna/RNA.v4.1.mix',

           # ---------- META ----------

          'data/meta/META.usage.txt',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v6.mix.csv',
          'data/meta/META.v1.tab.bed',

           # ---------- Variant ----------

          'data/var/DNA.usage.txt',
          'data/var/DNA.v3.mix.csv',
          'data/var/DNA.variant.bed',
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
