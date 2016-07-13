import os
import sys

def clean(x):
    x = x.replace("*", "");
    x = x.replace("/", "");
    x = x.replace("-", "");
    x = x.replace("+", "");
    x = x.replace("\n", "");    
    return (x)

#
# Usage: python varscan_to_Anaquin.py <varscan input file>
#

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print 'Usage: python varscan_to_Anaquin.py <varscan input file>'
    else:
        file = sys.argv[1]
    
        with open(file) as f:
            
            # Print the header
            print ('ChrID\tPosition\tRef\tAlt\tReadR\tReadV\tDepth\tQualR\tQualV\tPvalue')
                        
            i = 0
                        
            for line in f:
                i = i + 1
                
                # Skip the first line (header)
                if (i == 1):
                    continue
                
                toks = line.split('\t')
    
                # Does variant alelle included?
                if (len(toks) == 19):
                
                    # Eg: chrIS 
                    chrom = toks[0]
                
                    # Eg: 373703
                    pos  = toks[1]
            
                    # Eg: A
                    ref = toks[2]
                
                    # Eg: A
                    cons = toks[3]
        
                    # Eg: 102
                    reads1 = toks[4]
        
                    # Eg: 150
                    reads2 = toks[5]

                    # Variant frequency
                    varFreq = toks[6]
                
                    strs1   = toks[7]
                    strs2   = toks[8]
                    qual1   = toks[9]
                    qual2   = toks[10]
                    pval    = toks[11]
                    map1    = toks[12]
                    map2    = toks[13]
                    allele  = toks[18]
    
                    depth = reads1 + reads2
    
                    # Is this an insertion?
                    isInsert = '*/+' in cons
                
                    # Is this a deletion?
                    isDelete = '*/-' in cons
    
                    if (isInsert):
                        ref = clean(ref)
                        alt = clean(ref + allele)    
                    elif (isDelete):
                        ref = clean(ref + allele);
                        alt = clean(ref);
                    else:
                        ref = clean(ref)
                        alt = clean(allele)
                    
                    print (chrom + '\t' + pos + '\t' + ref + '\t' + alt + '\t' + reads1 + '\t' + reads2 + '\t' + depth + '\t' + qual1 + '\t' + qual2 + '\t' + pval)
