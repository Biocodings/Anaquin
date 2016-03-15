#!/usr/bin/python

#
# This script provides functionality for k-mer
#
#   index -i /Users/tedwong/Desktop/AFU007.v013.index /Users/tedwong/Desktop/AFU007.v013.fa
#   quant -i /Users/tedwong/Desktop/AFU007.v013.index -o /Users/tedwong/Desktop /Users/tedwong/Desktop/A.fq /Users/tedwong/Desktop/B.fq
#

import os
import sys

# The index used for consturcting k-mer names
kIndex = 1

class Sequin:
    
    # Eg: NG1_1_P1
    name = None
    
    # Length in bases
    length = None
    
    # The sequence of the sequin
    seq = None

class KMer:
    
    # Name of the K-mer
    name = None
    
    # The actual sequence
    seq = None

class PSL:
    
    # Name of the query
    qname = None
    
    # Name of the target
    tname = None
    
    # Start of the query
    qstart = None
    
    # End of the query
    qend = None

class Break:
    
    # Name of the chromosome
    name = None
    
    b1 = -1    
    b2 = -1

def readKMers(file):
    kmers = []
    with open(file) as f:
        while True:
            l = f.readline()
            if (not l):
                break
            if l[0] == '>':
                kmer = KMer()
                kmer.name = l[1:].strip() 
                kmer.seq  = '????'
                kmers.append(kmer)
    return kmers

def readPSL(file):    
    psls = {}
    with open(file) as f:
        while True:
            l = f.readline()
            if (not l):
                break
            if l[0] == '3':
                toks = l.split()
                
                p = PSL()

                p.qname  = toks[9]
                p.tname  = toks[13]                
                p.qstart = int(toks[15])
                p.qend   = int(toks[16])
                
                psls[p.qname] = p
    return psls

def buildKmers(seq, length):
    kmers = []
    
    i = 0
    j = length

    global kIndex

    while (j <= len(seq)):        
        kmer = KMer()
        kmer.name = 'K' +  str(kIndex)
        kmer.seq  = seq[i:j]                
        kmers.append(kmer)

        if len(kmer.seq) != length:
            raise Exception('len(kmer) != length')

        i = i + 1
        j = j + 1
        kIndex = kIndex + 1
        
    return kmers

def readStands(file, prefix):
    seqs = []
    with open(file) as f:
        while True:
            l = f.readline()
            if (not l):
                break
            
            toks = l.split()
                
            if len(toks) != 3:
                raise Exception('Invalid format: ' + file)
            
            if (toks[0].startswith(prefix)):
                seq = Sequin()
            
                seq.name   = toks[0]
                seq.seq    = toks[2]
                seq.length = int(toks[1])
        
                seqs.append(seq)

    return seqs    

def generateTempFA(kmers):
    out = open('target.fa', 'w')
    
    for kmer in kmers:
        out.write('>' + str(kmer.name) + '\n')
        out.write(kmer.seq + '\n')

    return 'target.fa'

if __name__ == '__main__':
    
    # Reference annotation (reference breakpoints)
    b_path = '/Users/tedwong/Sources/QA/data/FusQuin/AFU006.v032.bed'

    
    # Read in the normal genes
    nSeqs = readStands('/Users/tedwong/Sources/QA/data/FusQuin/AFU009.v032.csv', 'NG')

    # Read in the fusion genes
    fSeqs = readStands('/Users/tedwong/Sources/QA/data/FusQuin/AFU009.v032.csv', 'FG')

    if len(nSeqs) != len(fSeqs):
        raise Exception('len(nSeqs) != len(fSeqs)')

    #
    # 1. Read in the reference break points (relative to the chromosome)
    #

    breaks = {}

    with open(b_path) as f:
        while True:
            l = f.readline()
            if (not l):
                break
            toks = l.split()

            b = Break()
            
            b.name = toks[3]
            b.b1   = int(toks[1])
            b.b2   = int(toks[2])
    
            breaks[b.name] = b

    # The k-mer file we're generating
    out = open('AFU010.v032.csv', 'w')

    #
    # 2. Process normal genes
    #

    for nSeq in nSeqs:
        
        kmers = buildKmers(nSeq.seq, 31)
        query = generateTempFA(kmers)
        
        if nSeq.name != 'NG1_1_P2':
            continue
        
        # This is the only way we can find out whether the kmer is spanning...
        os.system('blat ' + 'CTR001.v021.fa' + ' ' + query + ' output.psl')

        # Read in the PSL alignment file
        psls = readPSL('output.psl')

        # This is the breakpoint we need
        b = breaks[nSeq.name]

        lastStr = None
        lastEnd = None
        
        # Is the k-mer in the spanning region?
        spanning = False

        # Number of spanning k-mers
        nSpans = 0

        for kmer in kmers:
            if kmer.name in psls:
                spanning = False
                psl = psls[kmer.name]                
                
                lastStr = psl.qstart
                lastEnd = psl.qend
                
                out.write(kmer.name + '\t' + 'NormN' + '\t' + str(psl.qstart) + '\t' + str(psl.qend) + '\n')
            else:
                
                isBackward = nSeq.name.endswith(('P2'))
                
                if not isBackward and lastEnd == b.b1:
                    spanning = True
                if isBackward and lastStr == b.b2:
                    spanning = True
                
                if spanning:
                    out.write(kmer.name + '\t' + 'NormS' + '\t' + '-' + '\t' + '-' + '\n')
                    nSpans = nSpans + 1
                else:
                    out.write(kmer.name + '\t' + 'NormN' + '\t' + '-' + '\t' + '-' + '\n')
        
        if nSpans == 0:
            raise Exception('No spanning k-mer found for ' + nSeq.name)
        
        
        
        
        
    print 'AFU010.v032.csv generated'        
    raise Exception('dsjkaskdsjkas')
        
        
    #
    # 3. Process fusion genes
    #

    for nSeq in nSeqs:
        kmers = buildKmers(nSeq.seq, 31)

        #
        # Now, we have the kmers and we need to know those kmers which span the fusion point.
        #


    out.close()
    



    raise Exception('kjsdjkadsjkds')
    
    
    

    # Where the sequence files are
    r_path = '/Users/tedwong/Sources/QA/data/FusQuin/AFU007.v013.fa'
    
    # Where the k-mers are
    k_path = '/Users/tedwong/Sources/QA/K'


    out = open('kmer.stats', 'w')

    #
    # Read in the k-mers, we'll need to do blat to generate the location in which the k-mers span across.
    # Assume that the k-mers have been generated by the C++ program.
    #
    #  Eg: blat /Users/tedwong/Sources/QA/data/FusQuin/AFU007.v013.fa NG1_12_P2 output.psl
    #
    # We'll find the spanning kmers for only and only this sequin
    #

    for seq in os.listdir(k_path):
        
        if (seq[0] != 'N'):
            continue

        if (seq != 'NG1_1_P1'):
            continue
            
        print ('Analyzing ' + seq)
            
        os.system('rm -rf output.psl')
        os.system('blat ' + 'CTR001.v021.fa' + ' ' + k_path + '/' + seq + ' output.psl')
        print (k_path + '/' + seq)

        # Now read the PSL file
        psls = readPSL('output.psl')

        #
        # Now we know where each of the K-mers are mapped to. We're only interested in the k-mers
        # crossing the fusion point.
        #
        
        # All the k-mers for the sequin
        kmers = readKMers(k_path + '/' + seq)

        for kmer in kmers:
            
            # Whether this k-mer can be found in the PSL file
            found  = False
            
            # What's the PSL alignment for this k-mer?
            for psl in psls:

                # Is this what we want?
                #if psl.qname == kmer and psl.tname == seq:
                if psl.qname == kmer.name:
                                        
                    # What's the known breakpoint for the sequin?
                    b = breaks[seq]
                    
                    # For simplicity, only the first breakpoint is considered... 
                    if psl.qstart <= b.b1 and psl.qend >= b.b1:
                        
                        out.write(kmer.name + '\t' + kmer.seq + '\t' + seq + '\t' + 'S2\n')
                        
                    else:
                        
                        # If we've found an alignment but it's not related to breakpoints
                        out.write(kmer.name + '\t' + kmer.seq + '\t' + seq + '\t' + 'S1\t' + str(psl.qstart) + '\t' + str(psl.qend) + 'ABBCD' + '\n')

                    found = True
                    break

            # This k-mer can't be found in the chromosome...
            if found == False:
                out.write(kmer.name + '\t' + kmer.seq + '\t' + seq + '\t' + 'S3' + '\n')
                #print ('Warning: PSL not found for: ' + kmer.name)
            #    pass

    out.close()







