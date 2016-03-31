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

          'src/r/plotMA.R',
          'src/r/plotSplice.R',
          'src/r/plotScatterPool.R',

          # ------- Fusion -------

          'src/r/plotROC_F.R',
          'src/r/plotExpress_F.R',

          'data/FusQuin/AFU004.v032.ref',
          'data/FusQuin/AFU005.v032.bed',
          'data/FusQuin/MFU007.v013.csv',

          # ------- Ladder -------

          'src/r/plotLadderAbund.R',

          'data/LadQuin/MLA014.v013.csv',
          'data/LadQuin/MLA016.v013.csv',
          'data/LadQuin/MLA020.v013.csv',
           
          # ---------- Transcriptome ----------

          'src/r/plotROC_T.R',
          'src/r/plotLODR_T.R',
          'src/r/plotScatter_T.R',

          'data/TransQuin/ATR001.v032.gtf',
          'data/TransQuin/MTR002.v013.csv',
          'data/TransQuin/MTR003.v013.csv',
          'data/TransQuin/MTR004.v013.csv',          
          'data/TransQuin/MTR005.v013.csv',

          # ---------- META ----------

          'data/MetaQuin/AME013.v032.fa',
          'data/MetaQuin/MME023.v013.csv',
          'data/MetaQuin/AME015.v032.bed',

          # ---------- Variant ----------

          'src/r/plotROC_V.R',
          'src/r/plotAlleleReads.R',
          'src/r/plotLODR_V.R',
          'src/r/plotAlleleAllele.R',
          'src/r/plotDensity.R',
          'src/r/plotSubsample.R',
	  'src/r/plotVAbundAbund.R',

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
