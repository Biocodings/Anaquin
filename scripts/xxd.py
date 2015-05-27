#
# This script converts resources into C++ headers such that the software can compile with.
#

#!/usr/bin/python

import os
import sys
import math
import subprocess

os.system('rm -rf src/data')

def xxd(src, dst):
    os.system('xxd -i ' + src + ' ' + ' ' + dst)

docs  = [ 'docs/manual.txt' ]
data  = [ 'data/silico.fa',
          'data/rna/standards.gtf',
          'data/rna/standards.bed',
          'data/rna/rna_standards.txt',
          'data/dna/variant.ChrT51.vcf',          
          'data/dna/DNA.ref.bed',
          'data/dna/DNA.tab.fa',
        ]
tests = [ ]

r = data + docs
for i in range(0,len(r)):
    file = os.path.basename(r[i])
    xxd(r[i], 'src/data/' + file)
