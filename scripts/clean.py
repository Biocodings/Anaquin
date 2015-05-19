#!/usr/bin/python

import os
import sys
import math
import subprocess

if __name__ == '__main__':
    w = open(sys.argv[2], 'w')

    with open(sys.argv[1]) as f:
        while True:
            l = f.readline()
            if (not l):
                break
            
            if l[0] == '@':
                w.write(l+'\n')
            else:
                if 'chrT' in l:
                    w.write(l+'\n')   
    w.close()