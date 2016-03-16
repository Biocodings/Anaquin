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
    
    # Is this an indel?
    isIndel = None
    
    # Position after reversing
    newPos = None

def readFA(file):
    parser = SeqIO.parse(open(file), 'fasta')
    x = {}

    for fasta in parser:
        name, seq = fasta.id, str(fasta.seq)
        
        f = FASTA()
        
        f.name = name
        f.seq  = seq

        x[f.name] = f
        
    return x

def readVCF(file):
    reader = vcf.Reader(open(file, 'r'))
    r = {}

    for record in reader:
        
        x = VCF()
        
        x.pos      = record.POS
        x.ref      = record.REF
        x.alt      = str(record.ALT[0])
        x.isSNP    = record.is_snp
        x.is_indel = record.is_indel
        
        
        x.refID = record.ID
        
        r['K' + str(x.pos)] = x

    return r

def uniToStr(i):
    return  'u\'\\u' + str(i)

#
# Usage: python seqtools reverse <genome.fa> <variants.vcf>
#
if __name__ == '__main__':

    # Eg: genome.fa
    rFile = 'test.fa'  #sys.argv[2]

    # Eg: variants.vcf
    vFile = 'test.vcf' #'/Users/tedwong/Sources/QA/data/VarQuin/AVA009.v032.vcf' #sys.argv[3]

    #
    # 1. Read in the reference genome
    #

    rData = readFA(rFile)
    
    if len(rData) > 1:
        raise Exception('More than a single chromosome is not supproted')

    #
    # 2. Read in the variant annotation
    #

    vData = readVCF(vFile)
    
    #
    # 3. Generate a alternative reference marker with GATK
    #

    rDict = os.path.splitext(rFile)[0] + '.dict'

    if False:
        print ('samtools faidx ' + rFile)

        print('Generating faidx index...\n')
        os.system('samtools faidx ' + rFile)
    
        print('Generating picard dictionary...\n')
        os.system('rm ' + rDict)
        os.system('java -jar /usr/lib/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=' + rFile + ' O=' + rDict)

        print('Generating FastaAlternateReferenceMaker...\n')
        os.system('java -jar /usr/lib/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ' + rFile + ' -o genome_with_variants.fa -V ' + vFile)
        print('\nFastaAlternateReferenceMaker genome (genome_with_variants.fa) generated\n')
    
    #
    # 4. Process the alternative reference (genome_with_variants.fa) by creating a mirror copy
    #
    
    aData = readFA('genome_with_variants.fa')
    
    # Reference sequence before reversing
    refSeq = aData.itervalues().next()
    
    # This is the sequence that will be substituted
    #altSeq = refSeq.seq
    
    # Length of the reference sequence
    n = len(refSeq.seq)
    
    #i = 2713
    #s = 'ABCD'
    
    #for j in range(10):        
    #    s = s + unichr(0x2713 + j)
            
    #print s
    #raise Exception('ddd')
    
    #
    # For every variant, replace all the entires with unicode so that once flipped, we know where to find them
    #
    
    for key in vData:
        
        if vData[key].isSNP:
            pos = vData[key].pos            
            print str(pos) + '\t' +  str(refSeq.seq[pos])

            #
            # What's the position of the SNP after reversing? It's n-(i+1)
            #

            vData[key].newPos = n-(pos+1)

        else:
            raise Exception('Indel not supported...')
    
    #
    # 5. Now, flip the genome...
    #
    
    flip = refSeq.seq[::-1]

    #
    # Generating the variant file
    #

    print ('\nWriting variant FASTA... Please check: [genome_with_vars.fa]')
    out = open('genome_with_vars.fa', 'w')
    out.write('>' + refSeq.name + '\n')
    out.write(refSeq.seq)
    out.close()

    #
    # Generating the flipped FASTA
    #
    
    print ('Writing flipped FASTA... Please check: [flipped_genome_with_vars.fa]')
    out = open('flipped_genome_with_vars.fa', 'w')
    out.write('>' + refSeq.name + '\n')
    out.write(flip)
    out.close()
    
    #
    # Generating the VCF file of variants relative to the flipped genome
    #
    
    print ('Writing VCF file... Please check: [flipped_vars.vcf]')    
    out = open('flipped_vars.vcf', 'w')
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
    out.write('##CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003\n')

    for key in vData:
        var = vData[key]        
        if var.isSNP:
            out.write(refSeq.name + '\t' + str(var.newPos) + '\t' + var.refID + '\t' + var.ref + '\t' + var.alt + '\n')
        else:
            raise Exception('Indel not supported...')            
    out.close()
    
    
    

    #samtools faidx test.fa
    #java -jar /usr/lib/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=CVA002.v015.fa O=CVA002.v015.dict
    #java -jar /usr/lib/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R CVA002.v015.fa -o ABCD.fa -V AVA009.v032.vcf







    
