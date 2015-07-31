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
           'data/fusion/fusion.bed',
           'data/fusion/normal.bed',                      

           # ------- Ladder -------

           'data/ladder/CON.v3.mix.csv',

           # ---------- Transcriptome ----------

          'data/trans/RNA.v1.gtf',
          'data/trans/RNA.v4.1.csv',

           # ---------- META ----------

          'data/meta/META.v1.tab.fa',
          'data/meta/META.v1.tab.fa',
          'data/meta/META.v6.mix.csv',
          'data/meta/META.v1.tab.bed',

           # ---------- Variant ----------

          'data/var/DNA.v3.mix.csv',
          'data/var/DNA.variant.bed',
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
