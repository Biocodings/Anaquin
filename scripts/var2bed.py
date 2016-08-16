import os
import sys

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        for line in f.readlines():
            line = line.strip()

            if line[0] != '#':
                toks  = line.split('\t')
                chrID = toks[0]
                pos = int(toks[1])
                
                start = pos - 500
                end   = pos + 500

                print (chrID + '\t' + str(start) + '\t' + str(end) + '\t' + 'D_' + str(start) + '_' + str(end))