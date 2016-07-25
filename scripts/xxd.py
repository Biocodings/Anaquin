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

          'scripts/reports.py',

          'src/r/plotFold.R',
          'src/r/plotScatter.R',
    	  'src/r/plotSensitivity.R',

          # ---------- Transcriptome ----------

          'src/r/plotTMinor.R',
          'src/r/plotTROC.R',
          'src/r/plotTLODR.R',          

          'data/RnaQuin/ARN020_v001.gtf',
          'data/RnaQuin/MRN027_v001.csv',
          'data/RnaQuin/MRN028_v001.csv',
          'data/RnaQuin/MRN029_v001.csv',          
          'data/RnaQuin/MRN030_v001.csv',

          # ---------- Variant ----------

   	      'src/r/plotVLOD.R',
          'src/r/plotVROC1.R',

          'data/VarQuin/sampled.bed',
          'data/VarQuin/AVA009_v001.vcf',
          'data/VarQuin/MVA011_v001.csv',
          'data/VarQuin/MVA012_v001.csv',
        ]
tests = [ ]

for file in os.listdir('data/manuals/'):
    path = os.path.join('data/manuals/', file)
    if os.path.isfile(path):
        os.system('expand -t 4 ' + path + ' > /tmp/tmp.txt')
        os.system('mv /tmp/tmp.txt ' + path)

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
                    
os.system('rm src/data/resources.o')