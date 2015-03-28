#!/usr/bin/python

import os
import sys
import math

os.system('rm -rf /home/tedwon/Projects/Temp')
os.system('mkdir /home/tedwon/Projects/Temp')
os.system('cd /home/tedwon/Projects/Temp')

os.system('cp /home/tedwon/Sources/QA/data/silico/silico.fa /home/tedwon/Projects/Temp')


os.system('wgsim -e0 -r0 -R0 silico.fa A.fastq B.fastq')



