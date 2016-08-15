import os
import sys
import numpy

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        c2l = {}
        
        for line in f.readlines():
            line = line.strip()

            if ',length=' in line:
                toks   = line.split(',length=')
                toks   = toks[1].split('>')
                length = toks[0]                
                toks   = line.split('ID=')
                toks   = toks[1].split(',length')
                chrID  = toks[0]                
                c2l[chrID] = int(length)
            
            if line[0] != '#':
                toks  = line.split('\t')
                chrID = toks[0]

                if 'chrev' in chrID:
                    pos  = toks[1]
                    pos  = c2l[chrID] - int(pos) + 1
                    
                    # Is the reference an indel?
                    isRefIndel = len(toks[3]) > 1
                    
                    # Is the alternative an indel?
                    isAltIndel = len(toks[4]) > 1
                    
                    if (isRefIndel):
                        pos = pos - len(toks[3])
                        
                    if (isAltIndel):
                        pos = pos - len(toks[4])

                    line = line.replace(chrID, chrID.replace('chrev', 'chr'))
                    line = line.replace(toks[1], str(pos))
                    
            print (line)    
            
