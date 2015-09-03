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

           'data/fusion/FUSMixture_3.0.csv',
           'data/fusion/FUSBreak_1.0.ref',
           'data/fusion/FUSFusionStandard_1.0.bed',
           'data/fusion/FUSFusionStandard_1.0.gtf',                      
           'data/fusion/FUSNormalStandard_1.0.bed',
           'data/fusion/FUSNormalStandard_1.0.gtf',

           # ------- Ladder -------

           'data/ladder/LadderMixture_3.0.csv',

           # ---------- Transcriptome ----------

          'data/trans/TransStandard_1.0.gtf',
          'data/trans/TransMixture_4.1.csv',

           # ---------- META ----------

          'data/meta/METAMixture_6.0.csv',
          'data/meta/METAStandard_1.0.bed',
          'data/meta/META_v1_tab.fa',

           # ---------- Variant ----------

          'data/var/AVA009.v032.vcf',
          'data/var/VARMixture_3.0.csv',
          'data/var/VARStandard_1.0.gtf',
          'data/var/VARVariant_1.0.bed',
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
