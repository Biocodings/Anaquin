#!/usr/bin/python

import os
import sys
import urllib

def run(cmd):
    os.system(cmd)

# Merge the inputs into a single SAM file
def merge(file, headers, humans, synthet):
    with open(file, 'w') as f:
        for l in headers:
            f.write(l)
        for l in synthet:
            f.write(l['line'])
        for l in humans:
            f.write(l['line'])

def create(file, headers, lines):
    with open(file, 'w') as f:
        for l in headers:
            f.write(l)
        for l in lines:
            f.write(l['line'])

def analyze(file):
    lengths = []  # Length of the reads
    headers = []  # Header lines    
    humans  = []  # Reads for the human genome
    synthet = []  # Reads for the synthetic chromosome

    # Number of bases sequenced for the human 
    l_humans = 0
    
    # Number of bases sequenced for the synthetic chromosome
    l_chrT = 0

    n = 0

    with open(file) as f:
        for l in f:
            n = n +1
            
            if (l[0] == '@'):
                headers.append(l)
            else:
                toks = l.split('\t')
                
                if (toks[2] == '*'):
                    continue
                else:
                    bases = len(toks[9])
                    
                    if (toks[2] == 'chrT'):
                        l_chrT = l_chrT + bases
                        synthet.append({ 'line': l, 'len': bases })
                    else:
                        l_humans = l_humans + bases
                        humans.append( { 'line': l, 'len': bases })

    return { 'headers': headers,
             'hg38'   : humans,
             'chrT'   : synthet }

def stats(file):
    r = analyze(file)

    # http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics
    hg38_len = 3209286105.0
    
    # grep -v chrT chrT.v0510.fa | wc -m
    chrT_len = 8457083.0

    # Estimate the coverage by:
    #
    #    number of mapped reads * read length
    #

    hg38_cov  = len(r['hg38']) * 125 / hg38_len
    chrT_cov  = len(r['chrT']) * 125 / chrT_len    
    subsample = hg38_cov / chrT_cov

    print('Number of reads for genome: '       + str(len(r['hg38'])))
    print('Number of reads for synthetic: '    + str(len(r['chrT'])))
    print('Estimated coverage for genome: '    + str(hg38_cov))
    print('Estimated coverage for synthetic: ' + str(chrT_cov))
    print('Need to subample in chrT: '         + str(subsample))    

# Normalize a SAM file to the human genome coverage
def normalize(input, output):
    r = analyze(input)

    # http://genomewiki.ucsc.edu/index.php/Hg38_7-way_Genome_size_statistics
    hg38_len = 3209286105.0
    
    # grep -v chrT chrT.v0510.fa | wc -m
    chrT_len = 8457083.0

    # Estimate the coverage by:
    #
    #    number of mapped reads * read length
    #

    hg38_cov  = len(r['hg38']) * 125 / hg38_len
    chrT_cov  = len(r['chrT']) * 125 / chrT_len    
    subsample = hg38_cov / chrT_cov

    print('Number of reads for genome: '       + str(len(r['hg38'])))
    print('Number of reads for synthetic: '    + str(len(r['chrT'])))
    print('Estimated coverage for genome: '    + str(hg38_cov))
    print('Estimated coverage for synthetic: ' + str(chrT_cov))
    print('Need to subample in chrT: '         + str(subsample))

    #
    # Now we have to estimated coverage for both genome and synthetic. We'll have to subsample the synthetic
    # reads to the coverage expected for the genome. Therefore, we'll randomly select some of the reads
    # aligned to chrT.
    #

    #
    # Create a SAM file for the chrT
    #
    
    create('__tmp__.sam', r['headers'], r['chrT'])
    run('samtools view -bS __tmp__.sam > __tmp__.bam')

    #
    # Ask samtools to do the subsampling
    #
    
    run('samtools view -s ' + str(subsample) + ' __tmp__.bam > __sampled__.sam')

    #
    # Create a new SAM file from the subsampled reads
    #
    
    create('__tmp1__.sam', r['headers'], r['hg38'])
    run('cat __sampled__.sam >> __tmp1__.sam')
    run('mv __tmp1__.sam ' + output)
    
if __name__ == '__main__':

    mode = sys.argv[1]

    # The raw SAM file  
    file = sys.argv[2]

    if (mode == '-n'):
        normalize(file, 'resampled.sam')
    elif (mode == '-s' or mode == 'stats'):
        stats(file)
        



