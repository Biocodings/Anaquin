#!/usr/bin/python

#
# This script provides subsampling to a SAM/BAM file for sequins
#

import os
import sys

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
    hg38 = []  # Reads for the human genome
    chrT = []  # Reads for the synthetic chromosome

    # Number of bases sequenced for the human 
    l_humans = 0
    
    # Number of bases sequenced for the synthetic chromosome
    l_chrT = 0

    with open(file) as f:
        for l in f:
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
                        chrT.append({ 'line': l, 'len': bases })
                    else:
                        l_humans = l_humans + bases
                        hg38.append( { 'line': l, 'len': bases })

    return { 'headers': headers, 'hg38' : hg38, 'chrT' : chrT }

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

    #
    # Estimate the coverage by:
    #
    #    number of mapped reads * read length / size of the chromosomes
    #

    hg38_cov = len(r['hg38']) * 125 / hg38_len
    chrT_cov = len(r['chrT']) * 125 / chrT_len
    
    # This is the ratio that we'll need to subsample chrT
    ratio = hg38_cov / chrT_cov

    print('Number of reads for genome: '       + str(len(r['hg38'])))
    print('Number of reads for synthetic: '    + str(len(r['chrT'])))
    print('Estimated coverage for genome: '    + str(hg38_cov))
    print('Estimated coverage for synthetic: ' + str(chrT_cov))
    print('Need to subample in chrT: '         + str(ratio))

    #
    # Let's subsample the synthetic reads to match the coverage expected for the genome. Therefore, we'll randomly sample without replacement for the reads
    # aligned to chrT.
    #

    #
    # Create a SAM file for the chrT
    #
    
    create('__tmp__.sam', r['headers'], r['chrT'])
    run('samtools view -bS __tmp__.sam > __tmp__.bam')

    #
    # Ask samtools to do the subsampling. Eg: samtools view -s <ratio> <BAM File> > <SAM File>
    #
    
    run('samtools view -s ' + str(ratio) + ' __tmp__.bam > __sampled__.sam')

    #
    # Merge the subsampled reads to create a new SAM/BAM file
    #

    create('__tmp1__.sam', r['headers'], r['hg38'])
    run('cat __sampled__.sam >> __tmp1__.sam')
    run('mv __tmp1__.sam ' + output)

#
# Eg: python subsample.py -n <SAM/BAM File>
#
if __name__ == '__main__':

    mode = sys.argv[1]

    # The raw SAM file  
    file = sys.argv[2]

    if (mode == '-n'):
        print('Normalizing ' + file)
        normalize(file, 'subsampled.sam')
    elif (mode == '-r'):
        print('Reporting ' + file)
        stats(file)



