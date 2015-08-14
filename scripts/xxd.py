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
           'data/fusion/FusionGenes.chrTv1.bed',
           'data/fusion/NormalParentGenes.chrTv1.bed',                      
           'data/fusion/FusionGenes.chrTv1.gtf',
           'data/fusion/NormalParentGenes.chrTv1.gtf',

           # ------- Ladder -------

           'data/ladder/Ladder_v3.csv',

           # ---------- Transcriptome ----------

          'data/trans/RNA_1.gtf',
          'data/trans/RNA_4_1.csv',

           # ---------- META ----------

          'data/meta/META_v6.csv',
          'data/meta/META_v1_tab.bed',
          'data/meta/META_v1_tab.fa',

           # ---------- Variant ----------

          'data/var/DNA_v3.csv',
          'data/var/DNA.variant.bed',
          'data/var/DNA.standards.chrT.gtf',
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
