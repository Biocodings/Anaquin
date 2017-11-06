#
# Python implementation for extracting unique k-mers for allele frequency ladder
#
#   python3 scripts/unique.py tests/data/A.V.23.fa A.V.8.vcf sequin_regions.combined.hg38.bed
#

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# BED file for the breakpoints
B1 = 'breaks.bed'

# BED file for k-mers
B2 = 'kmers.bed'

# K-mer length
k = 31

# Eg: CS_011, CS_012 ...
seqs = []

R = {}
V = {}

#
# Read BED file
#

BED = {}
with open(sys.argv[3]) as f:
    for line in f:
        toks = line.strip().split('\t')
        BED[toks[3]] = { 'start':int(toks[1]), 'end':int(toks[2]) }
assert(len(BED) > 0)

#
# Read FASTA file
#

recs = SeqIO.parse(sys.argv[1], 'fasta')
for r in recs:
    name = '_'.join(r.name.split('_')[0:2])
    
    if 'CS_' in name or 'CI_' in name:
        if '_R' in r.name:
            seqs.append(name)
            R[name] = str(r.seq)
        elif '_V' in r.name:
            seqs.append(name)
            V[name] = str(r.seq)

seqs = list(set(seqs))
seqs.sort()

assert(len(R) > 0)
assert(len(R) == len(V))

#
# Read VCF file
#

VCF = {}
with open(sys.argv[2]) as f:
    for line in f:
        if line[0] != '#':
            toks = line.strip().split('\t')
            
            ref = toks[3]
            alt = toks[4]
            
            isSNP = len(ref) == len(alt)
            isInd = not isSNP
            isIns = isInd and len(alt) > len(ref)
            isDel = isInd and not isIns
            isSub = isSNP            
            muts  = 'SNP' if isSNP else 'IND'
            type  = 'SUB' if isSub else 'INS' if isIns else 'DEL'

            assert(not toks[2] in VCF)
            VCF[toks[2]] = { 'name':toks[2], 'muts':muts, 'type':type, 'ref':ref , 'alt':alt, 'chr':toks[0], 'pos':int(toks[1]) }
assert(len(VCF) > 0)

#
# Let's find the breakpoints. It's simply the first unmatching base.
#

bs = {}
with open(B1, 'w') as f1:
    for seq in seqs:
        def __break__(s1, s2):
            i = 0
            while (s1[i] == s2[i]):
                i+=1
            return i

        # Find the breakpoint
        b = __break__(R[seq], V[seq])
        
        # What's the variant about?
        s = VCF[seq]['ref'] + '-' + VCF[seq]['alt']
        
        f1.write(str(seq) + '\t' + str(b-1) + '\t' + str(b) + '\t' + s + '\n')
        bs[seq] = b

#
# Now we can find out the starting and ending index for the unique k-mers
#

seqs_R = {}
seqs_V = {}

with open(B2, 'w') as f:
    for seq in seqs:
        R_ = R[seq]
        V_ = V[seq]
        I  = VCF[seq]

        # Breakpoint for the sequin
        b = bs[seq] + 1

        # Regular starting position
        rs = b-k
        
        # Unique k-mers for subtitution
        def SUB():
            return { 'startR':rs, 'endR':b+k-1, 'startV':rs, 'endV':b+k-1 }

        # Unique k-mers for insertion
        def INS():
            # Lenfth of allele
            l = len(I['ref'])
        
            # Regular ending
            re = b+k-(l-1)-1

            return { 'startR':rs, 'endR':b+k-1-(len(I['alt'])-1), 'startV':rs, 'endV':b+k-1 }

        # Unique k-mers for deletion
        def DEL():
            return { 'startR':rs, 'endR':b+len(I['ref'])-1+k-1, 'startV':rs, 'endV':b+len(I['alt'])-1+k-1 }

        if I['type'] == 'SUB':
            x = SUB()
        elif I['type'] == 'INS':
            x = INS()
        else:
            continue # FIX!
            x = DEL()

        sR = x['startR']
        eR = x['endR']
        sV = x['startV']
        eV = x['endV']
    
        f.write(seq + '_R\t' + str(sR) + '\t' + str(eR) + '\t' + seq + '_R\n')
        f.write(seq + '_V\t' + str(sV) + '\t' + str(eV) + '\t' + seq + '_V\n')
    
        #f.write(seq + '\t' + str(sR)   + '\t' + str(eR) + '\t' + seq + '_R\n')    # Debug
        #f.write(seq + '\t' + str(sV)   + '\t' + str(eV) + '\t' + seq + '_V\n')    # Debug
        #f.write(seq + '\t' + str(sR)   + '\t' + str(sR+k) + '\t' + seq + '_R1\n') # Debug
        #f.write(seq + '\t' + str(eR-k) + '\t' + str(eR) + '\t' + seq + '_R2\n')   # Debug
        #f.write(seq + '\t' + str(sV)   + '\t' + str(sV+k) + '\t' + seq + '_V1\n') # Debug
        #f.write(seq + '\t' + str(eV-k) + '\t' + str(eV) + '\t' + seq + '_V2\n')   # Debug

#
# Now we have the FASTA sequence file, and we have a BED file. Let's ask bedtools to do it for us.
#

os.system('bedtools getfasta -fi tests/data/A.V.23.fa -bed kmers.bed -fo kmers.fa')

#
# Let's check. All k-mers must be unique.
#

#
# Split a string into list of overlapping k-mers
#

def kmers(seq):
    kmers = []
    num_kmers = len(seq) - k + 1
    for i in range(num_kmers):
        kmer = seq[i:i+k]        
        kmers.append(kmer)        
    return kmers

x = {}
recs = SeqIO.parse('kmers.fa', 'fasta')
for r in recs:
    kms = kmers(str(r.seq))
    for km in kms:
        if km in x:
            print(Exception(r.name + ' has duplicate k-mer! (' + km + ')'))
        x[km] = 'True'
