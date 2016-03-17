#!/usr/bin/python

#
# This script provides functionality for k-mer
#
#   index -i /Users/tedwong/Desktop/AFU007.v013.index /Users/tedwong/Desktop/AFU007.v013.fa
#   quant -i /Users/tedwong/Desktop/AFU007.v013.index -o /Users/tedwong/Desktop /Users/tedwong/Desktop/A.fq /Users/tedwong/Desktop/B.fq
#

import os
import sys
import vcf
import time
from Bio import SeqIO

class FASTA:
    
    # Name of the chromosome
    name = None
    
    # Base sequenece
    seq = None

class VCF:
    
    refID = None
    
    # Eg: 7955610
    pos = None
    
    # Reference allele
    ref = None

    # Alternative alleles
    alt = None
    
    # Is this a SNP?
    isSNP = None
    
    # Is this an insertion?
    isInsert = False
    
    # Is this a deletion?
    isDelete = False
    
    # Position after flipping
    newPos = None
    
    # This is used for tracking the original position in the input file
    origPos = None

def readFA(file, chromID):
    parser = SeqIO.parse(open(file), 'fasta')
    x = {}

    for fasta in parser:
        name, seq = fasta.id, str(fasta.seq)
        
        if name == chromID:
            f = FASTA()
        
            f.name = name
            f.seq  = seq

            x[f.name] = f
        
    return x

def readVCF(file):
    reader = vcf.Reader(open(file, 'r'))
    r = []

    for record in reader:        
        x = VCF()
        
        x.pos      = record.POS
        x.ref      = record.REF
        x.refID    = record.ID
        x.origPos  = record.POS
        x.alt      = str(record.ALT[0])
        x.isSNP    = record.is_snp
        x.is_indel = record.is_indel
        
        r.append(x)

    return r

def uniToStr(i):
    return  'u\'\\u' + str(i)

#
# Custom implementation of FastaAlternateReferenceMaker in GATK
#
def FastaAlternateReferenceMaker(seq, vcfs):

    #
    # 1. Sort and classify the VCFs by position and types
    # 2. Go through the indels by reverse position and update structures (anything after the position)
    # 3. Use the updated information for SNP
    #
    
    print ('Number of variants: ' + str(len(vcfs)))
    vcfs.sort(key=lambda x: x.pos, reverse=True)
    print ('Variants sorted!')

    i = 0

    #
    # Process indels (they're processed first because they'd change the location of all variants)
    #

    for vcf in vcfs:
        print str(i) + '/' + str(len(vcfs))
        i = i + 1
        
        if vcf.isSNP:
            continue
        
        # The number depends on the variant context
        n = len(vcf.alt) - len(vcf.ref)

        if vcf.isInsert:
            seq = seq[0:vcf.pos] + vcf.alt + seq[vcf.pos+1:]
        else:
            seq = seq[0:vcf.pos] + vcf.alt + seq[vcf.pos + abs(n) + 1:]
    
        # Since we always iterate from backward, only the variants after it needs to be updated
        for t in vcfs:            
            # Does the variant occur after?
            if t.pos > vcf.pos:
                t.pos = t.pos + n

    #
    # Process SNPs now... This is easy because no position is changed
    #

    i = 0

    for vcf in vcfs:
        print str(i) + '/' + str(len(vcfs))
        i = i + 1

        if vcf.isSNP:
            seq = seq[0:vcf.pos] + vcf.alt + seq[vcf.pos+1:]

    return seq, vcfs

def createTestData():
    #
    # Say we have the following sequence: 0123456789ABCDEFGHIJKLMNOPQ and:
    #
    #    - SNP:       pos 5, from '5' to '#'
    #    - Insertion: pos 9, from '9' to @@@@
    #    - Deletion:  pos 13, DEFGH to &
    #    - SNP:       pos 26, from 'Q' to '!'
    #
    #  01234#678@@@@ABC&IJKLMNOP!
    #  !PONMLKJI&CBA@@@@876#43210
    #

    v1 = VCF()
    v1.pos = 5
    v1.ref = '5'
    v1.alt = '#'
    v1.isSNP = True

    v2 = VCF()
    v2.pos = 9
    v2.ref = '9'
    v2.alt = '@@@@'
    v2.isInsert = True

    v3 = VCF()
    v3.pos = 13
    v3.ref = 'DEFGH'
    v3.alt = '&'
    v3.isDelete = True
    
    v4 = VCF()
    v4.pos = 26
    v4.ref = 'Q'
    v4.alt = '!'
    v4.isSNP = True

    return '0123456789ABCDEFGHIJKLMNOPQ', [v1, v2, v3, v4]

def testSubtitute():
    
    seqs, vcfs = createTestData()    
    seqs, vcfs = FastaAlternateReferenceMaker(seqs, vcfs)

    # 01234#678@@@@ABC&IJKLMNOP!
    print seqs

def generateVCF(file, vcfs):
    out = open(file, 'w')
    out.write('##fileformat=VCFv4.0\n')
    out.write('##fileDate=20090805\n')
    out.write('##source=myImputationProgramV3.1\n')
    out.write('##reference=1000GenomesPilot-NCBI36\n')
    out.write('##phasing=partial\n')
    out.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
    out.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
    out.write('##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">\n')
    out.write('##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n')
    out.write('##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">\n')
    out.write('##CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT\n')

    for vcf in vcfs:
        out.write(chromID + '\t' + str(vcf.newPos) + '\t' + str(vcf.origPos) + '\t' + vcf.ref + '\t' + vcf.alt + '\n')

    out.close()

#
# Usage: python seqtools reverse <genome.fa> <variants.vcf>
#
# Eg: python seqtools reverse 
#
if __name__ == '__main__':

    # Eg: chr21
    chromID = sys.argv[1]

    # Eg: genome.fa
    rFile = sys.argv[2]

    # Eg: variants.vcf
    vFile = sys.argv[3]

    #
    # 1. Read in the variant annotation
    #
    start = time.time() 
    vData = readVCF(vFile)
    end   = time.time()
    
    print str(end - start) + 's taken for reading variants'

    #
    # 3. Generate a alternative reference marker with GATK
    #

    rDict = os.path.splitext(rFile)[0] + '.dict'

    if False:

        print('Generating faidx index...\n')
        
        # Eg: samtools faidx Hg38.clean.fa
        os.system('samtools faidx ' + rFile)
    
        print('Generating picard dictionary...\n')
        os.system('rm ' + rDict)
        
        # Eg: java -jar /usr/lib/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=Hg38.clean.fa O=Hg38.clean.dict
        os.system('java -jar /usr/lib/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=' + rFile + ' O=' + rDict)

        print('Generating FastaAlternateReferenceMaker...\n')

        #
        # Due to performance, this should be done by the user outside of the script
        #

        # Eg: java -jar /usr/lib/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R Hg38.clean.fa -V test.vcf
        os.system('java -jar /usr/lib/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + rFile + ' -o genome_with_variants.fa -V ' + vFile)

    #
    # 2. Read in the reference genome. This would represent the perfect genome.
    #

    start = time.time() 
    rData = readFA(rFile, chromID)
    end   = time.time()
    
    if len(rData) == 0:
        raise Exception('No chromosome is found: ' + chromID)
    if len(rData) > 1:
        raise Exception('More than a single chromosome is not supproted')

    print str(end - start) + 's taken for reading the genome'

    #
    # 3. Create the variant genome. This would be the genome representing high-confident variants. We should do it
    #    per chromosome.
    #

    rSeq = rData[chromID].seq
    vDat = vData

    #rSeq, vDat = createTestData()

    print('Running Mercer FastaAlternateReferenceMaker')
    vSeq, vcfs = FastaAlternateReferenceMaker(rSeq, vDat)
    print('Completed Mercer FastaAlternateReferenceMaker')
    
    if len(vcfs) != len(vDat):
        print len(vcfs)
        print len(vData)
        raise Exception('len(vcfs) != len(vData)')
    
    # Length of the reference sequence
    n = len(vSeq)
    
    print ('Number of bases: ' + str(n))
    
    #
    # For every variant, calculate the new position after flipped. This is possible because we know the number of charact
    #
    
    for vcf in vcfs:
        pos = vcf.pos            
        
        if vcf.isSNP:
            # What's the position of the SNP after reversing? It's n-(i+1)
            vcf.newPos = n-(pos+1)

        elif vcf.isInsert:
            
            # Number of characters inserted
            offset = len(vcf.alt) - len(vcf.ref)

            vcf.newPos = n-(pos+1)-offset

        else:
            vcf.newPos = n-(pos+1)
            
    #
    # 4. Now, flip the genome...
    #

    print 'Flipping the genome...'
    flip = vSeq[::-1]
    print 'Flipping completed'

    #
    # Generating the variant file
    #

    print ('\nWriting variant FASTA... Please check: [genome_with_vars.fa]')
    out = open('genome_with_vars.fa', 'w')
    out.write('>' + chromID + '\n')
    out.write(vSeq)
    out.close()

    #
    # Generating the flipped FASTA
    #
    
    print ('Writing flipped FASTA... Please check: [flipped_genome_with_vars.fa]')
    out = open('flipped_genome_with_vars.fa', 'w')
    out.write('>' + chromID + '\n')
    out.write(flip)
    out.close()
    
    #
    # Generating the VCF file of variants relative to the flipped genome
    #

    print ('Writing VCF file... Please check: [flipped_vars.vcf]')    
    generateVCF('flipped_vars.vcf', vcfs)