import os
import sys

#
# Convert a synthetic VCF file, add flanking regions to each of the variants and generate a new BED file.
#

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        for line in f.readlines():
            line = line.strip()

            if line[0] != '#':
                toks  = line.split('\t')
                chrID = toks[0]
                seqID = toks[2]
                pos = int(toks[1])
                
                start = pos - 500
                end   = pos + 500

                print (chrID + '\t' + str(start) + '\t' + str(end) + '\t' + seqID)