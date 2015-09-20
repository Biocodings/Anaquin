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
          'scripts/viewer.py',

           # ------- Fusion -------

           'data/fusion/AFU004.v032.ref',
           'data/fusion/AFU005.v032.bed',
           'data/fusion/MFU007.v013.csv',

           # ------- Ladder -------

           'data/ladder/LadderMixture_3.0.csv',

           # ---------- Transcriptome ----------

          'data/trans/ATR001.v032.gtf',
          'data/trans/MTR002.v013.csv',
          'data/trans/MTR003.v013.csv',
          'data/trans/MTR004.v013.csv',          

           # ---------- META ----------

          'data/meta/META_v1_tab.fa',
          'data/meta/MME023.v013.csv',
          'data/meta/METAStandard_1.0.bed',

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
