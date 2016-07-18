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
          'data/manuals/RnaAlign.txt',
          'data/manuals/RnaAssembly.txt',
          'data/manuals/RnaExpression.txt',
          'data/manuals/RnaFoldChange.txt',
          'data/manuals/RnaSubsample.txt',

          'data/manuals/VarAlign.txt',
          'data/manuals/VarSubsample.txt',
          'data/manuals/VarDiscover.txt',

          'scripts/viewer.py',
          'scripts/reports.py',
          'scripts/sleuth.R',

          'src/r/plotFold.R',
          'src/r/plotScatter.R',
    	  'src/r/plotSensitivity.R',

          # ------- Fusion -------

          'src/r/plotFROC.R',
          'src/r/plotFNormal.R',
	      'src/r/plotFFusion.R',
	      'src/r/plotFFold.R',

          'data/FusQuin/AFU004.v032.bed',
          'data/FusQuin/AFU005.v032.bed',
          'data/FusQuin/MFU007.v013.csv',

          # ------- Ladder -------

          'src/r/plotLNorm.R',

          'data/LadQuin/MLA014.v013.csv',
          'data/LadQuin/MLA016.v013.csv',
          'data/LadQuin/MLA020.v013.csv',
           
          # ---------- Transcriptome ----------

          'src/r/plotTMA.R',
          'src/r/plotTMinor.R',
          'src/r/plotTROC.R',
          'src/r/plotTLODR.R',          

          'data/RnaQuin/ARN020_v001.gtf',
          'data/RnaQuin/MRN027_v001.csv',
          'data/RnaQuin/MRN028_v001.csv',
          'data/RnaQuin/MRN029_v001.csv',          
          'data/RnaQuin/MRN030_v001.csv',

          # ---------- META ----------

          'src/r/plotMFold.R',
          'src/r/plotMKMer.R',
          'src/r/plotMReads.R',
          'src/r/plotMKAbund.R',
          'src/r/plotMAssembly.R',

          'data/MetaQuin/AME013.v032.fa',
          'data/MetaQuin/MME023.v013.csv',
          'data/MetaQuin/AME015.v032.bed',

          # ---------- Variant ----------

   	      'src/r/plotVLOD.R',
          'src/r/plotVROC1.R',
          'src/r/plotVROC2.R',          
          'src/r/plotVDensity.R',

          'data/VarQuin/sampled.bed',
          'data/VarQuin/AVA009_v001.vcf',
          'data/VarQuin/MVA011_v001.csv',
          'data/VarQuin/MVA012_v001.csv',
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