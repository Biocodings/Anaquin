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
          'scripts/viewer.py',

          'src/r/plotLODR.R',
          'src/r/plotScatter.R',

           # ------- Fusion -------

          'data/fusion/AFU004.v032.ref',
          'data/fusion/AFU005.v032.bed',
          'data/fusion/MFU007.v013.csv',

           # ------- Ladder -------

          'data/ladder/MLA014.v013.csv',
          'data/ladder/MLA016.v013.csv',
          'data/ladder/MLA020.v013.csv',
           
           # ---------- Transcriptome ----------

          'data/trans/ATR001.v032.gtf',
          'data/trans/MTR002.v013.csv',
          'data/trans/MTR003.v013.csv',
          'data/trans/MTR004.v013.csv',          
          'data/trans/MTR005.v013.csv',

           # ---------- META ----------

          'data/meta/AME013.v032.fa',
          'data/meta/MME023.v013.csv',
          'data/meta/AME015.v032.bed',

           # ---------- Variant ----------

          'data/var/AVA008.v032.bed',
          'data/var/AVA009.v032.vcf',
          'data/var/MVA011.v013.csv',
          'data/var/MVA012.v013.csv',
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
