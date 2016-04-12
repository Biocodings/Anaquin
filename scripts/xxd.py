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
          'scripts/reports.py',
          'scripts/sleuth.R',

          # ------- Fusion -------

          'src/r/plotFROC.R',
          'src/r/plotFExpress.R',

          'data/FusQuin/AFU004.v032.ref',
          'data/FusQuin/AFU005.v032.bed',
          'data/FusQuin/MFU007.v013.csv',

          # ------- Ladder -------

          'src/r/plotLExpress.R',

          'data/LadQuin/MLA014.v013.csv',
          'data/LadQuin/MLA016.v013.csv',
          'data/LadQuin/MLA020.v013.csv',
           
          # ---------- Transcriptome ----------

          'src/r/plotTMA.R',
	      'src/r/plotTFold.R',
	      'src/r/plotTMajor.R',
          'src/r/plotTROC.R',
          'src/r/plotTLODR.R',          
	      'src/r/plotTMultiple.R',
          'src/r/plotTExpress.R',

          'data/TransQuin/ATR001.v032.gtf',
          'data/TransQuin/MTR002.v013.csv',
          'data/TransQuin/MTR003.v013.csv',
          'data/TransQuin/MTR004.v013.csv',          
          'data/TransQuin/MTR005.v013.csv',

          # ---------- META ----------

          'src/r/plotMFold.R',

          'data/MetaQuin/AME013.v032.fa',
          'data/MetaQuin/MME023.v013.csv',
          'data/MetaQuin/AME015.v032.bed',

          # ---------- Variant ----------

          'src/r/plotVROC.R',
          'src/r/plotVAlleleReads.R',
          'src/r/plotVProb.R',
          'src/r/plotVAllele.R',
          'src/r/plotVDensity.R',
          'src/r/plotVSubsample.R',
	      'src/r/plotVExpress.R',

          'data/VarQuin/AVA017.v032.bed',
          'data/VarQuin/AVA009.v032.vcf',
          'data/VarQuin/MVA011.v013.csv',
          'data/VarQuin/MVA012.v013.csv',
        ]
tests = [ ]

r = data
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/resources/' + file)
    
    for file in os.listdir('src/resources/'):
        path = os.path.join('src/resources/', file)
        if os.path.isfile(path):
            with open(path, 'r') as f:
                data = f.read()
                if ('0x0a, 0x0a\n' in data):
                    raise Exception('Error: ' + path)
