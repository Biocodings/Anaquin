#!/usr/bin/python

import os
import random

def readMixA():
    
    file  = '/home/tedwon/Sources/QA/data/ladder/MLA015.v013.csv'
    file  = '/Users/tedwong/Sources/QA/data/ladder/MLA015.v013.csv'

    mixes = {}
    
    # The scaling factor
    fold = 787
    
    with open(file, 'r') as f:
        for line in f:            
            # Skip the first line
            if line[0] != 'I':
                toks = line.split('\t')
                
                name =  toks[0]
                conc =  float(toks[1])
                conc *= fold
                
                if round(conc) < 1:
                    raise Exception('conc < 1')

                mixes[name] = int(round(conc))

    return mixes
    
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches

def readSequence(file):
    seq = ''
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] != '>':
                seq = seq + line
    return seq

#
# The fourth generation of the simulator. This one was the version used for /tim/mer/CNV/Ira.
#

def simulate_4():

    print ('Simulation 4')

    # Read in for chromosome 21
    chr21 = readSequence('chr21.fa')
    
    print ('Number of bases in the original: ' + str(len(chr21)))

    # Read in for D4Z4.fa
    d4z4 = readSequence('Ira/D4Z4.fa')

    # Read in for chrM.fa
    chrM = readSequence('Ira/chrM.fa')

    # Read in for rDNA.fa
    rDNA = readSequence('Ira/rDNA.fa')

    seqs = []
    seqs.append(d4z4)
    seqs.append(chrM)
    seqs.append(rDNA)

    names = []
    names.append('d4z4')
    names.append('chrM')
    names.append('rDNA')

    chrS = open('Simulate4/chrS.fa',  'w')
    bedS = open('Simulate4/chrS.bed', 'w')

    fold = 20
    
    # We'll need another copy for chr21 to restore from it...
    chr21_ = chr21

    #
    # Generate a reference genome by adding a copy of each of standard once and only once
    #

    # Find a locus randomly for d4z4
    d4z4_i = random.randrange(22000000, len(chr21))
        
    # Find a locus randomly for chrM
    chrM_i = random.randrange(22000000, len(chr21))
        
    # Find a locus randomly for rDNA
    rDNA_i = random.randrange(22000000, len(chr21))
        
    chr21 = chr21[:d4z4_i] + d4z4 + chr21[d4z4_i:]
    chr21 = chr21[:chrM_i] + chrM + chr21[chrM_i:]
    chr21 = chr21[:rDNA_i] + rDNA + chr21[rDNA_i:]

    # What's the location for d3z4?
    d3z4_i = chr21.find(d4z4)
    
    # What's the location for d3z4?
    chrM_i = chr21.find(chrM)

    # What's the location for d3z4?
    rDNA_i = chr21.find(rDNA)
    
    if (d3z4_i == -1 or chrM_i == -1 or rDNA_i == -1):
        raise Exception('?????????')

    bedR = open('Simulate4/chrR.bed', 'w')
    bedR.write('chrR' + '\t' + str(d3z4_i) + '\t' + str(d3z4_i + len(d4z4))  + '\t' + 'd3z4' + '\n')    
    bedR.write('chrR' + '\t' + str(chrM_i) + '\t' + str(chrM_i + len(chrM))  + '\t' + 'chrM' + '\n')    
    bedR.write('chrR' + '\t' + str(rDNA_i) + '\t' + str(rDNA_i + len(rDNA))  + '\t' + 'rDNA' + '\n')    
    bedR.close()

    chrR = open('Simulate4/chrR.fa',  'w')
    chrR.write('>chrR\n')    
    chrR.write(chr21 + '\n')
    chrR.close()

    #
    # Generate random copies all over the genome...
    #

    chr21 = chr21_

    # Randomly insert all standards into the genome (repeat for n times)
    for n in range(fold):
        print ('Generating for fold: ' + str(n+1))
        
        # Find a locus randomly for d4z4
        d4z4_i = random.randrange(22000000, len(chr21))
        
        # Find a locus randomly for chrM
        chrM_i = random.randrange(22000000, len(chr21))
        
        # Find a locus randomly for rDNA
        rDNA_i = random.randrange(22000000, len(chr21))
        
        chr21 = chr21[:d4z4_i] + d4z4 + chr21[d4z4_i:]
        chr21 = chr21[:chrM_i] + chrM + chr21[chrM_i:]
        chr21 = chr21[:rDNA_i] + rDNA + chr21[rDNA_i:]
        
    print ('Finding the insertion points')

    # Trying to find all the insertion points... (the chromosome is not modified here)
    for i in range(3):
        matches = list(find_all(chr21, seqs[i]))

        if len(matches) == 0:
            raise Exception('Unknown???')            

        for j in range(len(matches)):

            # This is the starting position of the sequence
            start = matches[j]

            # This is the ending position of the sequence
            end = matches[j] + len(seqs[i])

            bedS.write('chrR' + '\t' + str(start) + '\t' + str(end)  + '\t' + names[i] + '\n')

    print ('Final number of bases in CNV: ' + str(len(chr21)))

    #
    # Generating a chromosome for simualtion
    #

    chrS.write('>chrR\n')
    chrS.write(chr21 + '\n')
    chrS.close()

    bedS.close() 

def generateReference(chr21, names, seqs, faFile, bedFile):

    for i in range(len(seqs)):
        
        # Find a locus randomly for the sequence
        #r =  random.randrange(22000000, len(chr21))
        r = 22000000 + (5000000 * (i+1))
         
        # Insert the sequence...
        chr21 = chr21[:r] + seqs[i] + chr21[r:]

    bedR = open(bedFile, 'w')

    for i in range(len(seqs)):

        # What's the location?
        r = chr21.find(seqs[i])

        if (r == -1):
            raise Exception('Failed in generateReference()')

        bedR.write('chrR' + '\t' + str(r) + '\t' + str(r + len(seqs[i]))  + '\t' + names[i] + '\n')

    bedR.close()

    chrR = open(faFile,  'w')
    chrR.write('>chrR\n')    
    chrR.write(chr21 + '\n')
    chrR.close()
    
    print ('Reference generated')

# Generate simulation for mixture A
def generateSimulationForA(chr21, names, seqs, faFile, bedFile):

    mixes = readMixA()

    chrS = open(faFile,  'w')
    bedS = open(bedFile, 'w')

    #
    # Generate random copies all over the genome...
    #

    k = 0

    for i in range(len(names)):
        print 'Generating for ' + names[i]
            
        fold = mixes[names[i]]

        if (fold >= 500):
            print 'Skip ' + names[i]
            continue

        for j in range(fold):

            # Find a locus randomly for the sequence
            r = 32000000 + (70000 * (k+1))
            #r =  random.randrange(22000000, len(chr21))

            # Insert the sequence...
            chr21 = chr21[:r] + seqs[i] + chr21[r:]
            
            # Make sure we don't overlap the sequences...
            k = k + 1

    #
    # Trying to find all insertion points and generate an annotation file
    #
            
    print ('Finding the insertion points')

    for i in range(len(seqs)):
        matches = list(find_all(chr21, seqs[i]))

        if len(matches) == 0:
            print ('Unknown to find: ' + names[i] + ' . Skipped.')
            continue

        for j in range(len(matches)):
            # This is the starting position of the sequence
            start = matches[j]

            # This is the ending position of the sequence
            end = matches[j] + len(seqs[i])

            bedS.write('chrR' + '\t' + str(start) + '\t' + str(end)  + '\t' + names[i] + '\n')    

    bedS.close()

    #
    # Generating a chromosome for simulation
    #

    chrS.write('>chrR\n')
    chrS.write(chr21 + '\n')
    chrS.close()

def generateSimulation(chr21, fold, names, seqs, faFile, bedFile):

    chrS = open(faFile,  'w')
    bedS = open(bedFile, 'w')

    #
    # Generate random copies all over the genome (repeat for n time)...
    #

    k = 0

    for n in range(fold):
        print ('Generating for fold: ' + str(n+1))
        
        for i in range(len(seqs)):
             
            # Find a locus randomly for the sequence
            r = 32000000 + (50000 * (k+1))
            #r =  random.randrange(22000000, len(chr21))
        
            # Insert the sequence...
            chr21 = chr21[:r] + seqs[i] + chr21[r:]
            
            # Make sure we don't overlap the sequences...
            k = k + 1
        
    print ('Finding the insertion points')

    # Trying to find all the insertion points... (the chromosome is not modified here)
    for i in range(len(seqs)):
        matches = list(find_all(chr21, seqs[i]))

        if len(matches) != fold:
            raise Exception('!!!!!!!')            

        for j in range(len(matches)):

            # This is the starting position of the sequence
            start = matches[j]

            # This is the ending position of the sequence
            end = matches[j] + len(seqs[i])

            bedS.write('chrR' + '\t' + str(start) + '\t' + str(end)  + '\t' + names[i] + '\n')    

    bedS.close()

    #
    # Generating a chromosome for simualtion
    #

    chrS.write('>chrR\n')
    chrS.write(chr21 + '\n')
    chrS.close()

#
# This version generates standards for the entire genome according to the fold changes.
#

def simulate_5():
    
    print ('Running simulation 5')
    
    # Read in for chromosome 21
    chr21 = readSequence('chr21.fa')

    print ('Number of bases in the original: ' + str(len(chr21)))

    #
    # Read in all the standards
    #

    stands = None
    with open('repeat_inserts_for_ted.csv', 'r') as f:
        stands = f.readlines()
        
    names = []
    seqs  = []
    
    for i in range(len(stands)):
        toks = stands[i].split(',')
        names.append(toks[0].strip())
        seqs.append(toks[1].strip())

    generateReference(chr21, names, seqs, 'Simulate5/chrR.fa', 'Simulate5/chrR.bed')
    generateSimulationForA(chr21, names, seqs, 'Simulate5/chrS.fa', 'Simulate5/chrS.bed')

#
# This version generates the standards for the entire genome (flat mix) and is successfully.
#

def simulate_3():
    
    # Read in for chromosome 21
    chr21 = readSequence('chr21.fa')

    print ('Number of bases in the original: ' + str(len(chr21)))

    #
    # Read in all the standards
    #

    stands = None
    with open('repeat_inserts_for_ted.csv', 'r') as f:
        stands = f.readlines()
        
    names = []
    seqs  = []
    
    for i in range(len(stands)):
        toks = stands[i].split(',')
        names.append(toks[0].strip())
        seqs.append(toks[1].strip())

    fold = 20

    generateReference(chr21, names, seqs, 'Simulate3/chrR.fa', 'Simulate3/chrR.bed')
    generateSimulation(chr21, fold, names, seqs, 'Simulate3/chrS.fa', 'Simulate3/chrS.bed')

def whereAreTheReads():
    
    # Where we want to filter our simulated reads
    i = 30412205

    # Where we want to filter our simulated reads
    j = 30455201

    # Name of the alignment file
    file = 'Simulate4/aligned_sorted.sam'

    with open(file, 'r') as f:        
        for line in f:
            line = line.strip()            
            toks = line.split('\t')             
            
            # Where does this read mapped to?
            mapped = toks[7]
            
            toks = toks[0]
            
            # Eg: chrF_2_536_1:0:0_1:0:0_16f881c
            toks = toks.split('_')
            
            start = int(toks[1])
            end   = int(toks[2])
            
            # This is the read we're interested, where does that mapped to?
            if (start >= i and end <= j):
                print mapped
                
            raise Exception('dddd')

def simulate():
    
    # Where we'll insert standard into (usually chr21)
    refFile = 'CVA002.v015.fa'
    
    # How many copies for each of the standard?
    fold = 10
    
    # Set this if only want a single standard be generated
    onlyFirst = True
    
    # Set this if we want to copy standards from the genome
    copyFrom = False
    
    #
    # 1. Read the sequence file for the original genome
    #
    
    chr21 = ''
    
    with open(refFile, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] != '>':
                chr21 = chr21 + line

    print ('Number of bases in the original: ' + str(len(chr21)))

    chrR = open('Simulate/chrR.fa',  'w')
    chrS = open('Simulate/chrS.fa',  'w')    
    bedR = open('Simulate/bedR.bed', 'w')
    bedS = open('Simulate/bedS.bed', 'w')

    #
    # 2. Generate a modified genome (either copying or inserting)
    #

    # We'll modify this reference genome if necessary
    ref = chr21

    # The copies that we'll insert
    copies = None

    j = None

    # Where to start?
    i = 10673

    if copyFrom:
        # How many characters to copy from?
        n = 600
    
        j = i + n
    
        # This is our sequence...
        copies = ref[i:j]
        
        # Starting for the copied sequence
        start = i
    
        # Ending fo the copied sequence
        end = j
    
        bedR.write('chrR\t' + str(i) + '\t' + str(j) + '\t' + 'Copied\n')    
        bedS.write('chrR\t' + str(i) + '\t' + str(j) + '\t' + 'Copied\n')        

        print ('Number of bases in CNV: ' + str(len(copies)))
        print ('The initial sequences in the copy is: ' + copies[:10])

    else:
        stands = None
        with open('repeat_inserts_for_ted.csv', 'r') as f:
            stands = f.readlines()

            names = []
            seqs  = []
    
        for i in range(len(stands)):
            toks = stands[i].split(',')
            names.append(toks[0].strip())
            seqs.append(toks[1].strip())

            # Randomly insert all standards into the genome (repeat for n times)
            for n in range(fold):
                print ('Generating for fold: ' + str(n+1))
                for k in range(len(names)):
                    #i = random.randrange(1, len(chr21))
                    i = (2*j) + (n * 10000) # Separated by 10000...
                    chr21 = chr21[:i] + '@' + str(k) + '#' + seqs[k] + chr21[i:]
                    if onlyFirst:
                        break

            # Very important, only generate a reference for the first fold (all other folds simulate CNV)
            if n == 1:            
                chr21_ = chr21

    #
    # Now we've got the reference genome, we'll insert copies into it.
    #

    # This is the genome for simulation
    sim = ref

    for f in range(fold):
        print ('Generating for fold: ' + str(f+1))
        for k in range(len(names)):
            ki = (2*j) + (k * 10000) # Separated by 10000...
            sim = sim[:ki] + '@' + str(k) + '#' + copies + sim[ki:]
            if onlyFirst:
                break

    #
    # Find all the insertion points
    #

    #i = 1
    #for k in range(len(sim)):
    #    if sim[k] == '@':
    #        bedS.write('chrR\t' + str(k+1) + '\t' + str(k + n + 1) + '\t' + 'chrF_' + str(i) + '\n')
    #        i = i + 1

    # Trying to find all the insertion points... (the chromosome is not modified here)
    for i in range(len(names)):
        matches = list(find_all(chr21, '@' + str(i) + '#'))

        if len(matches) == 0 and onlyFirst == False:
            raise Exception('Unknown???')            

        for j in range(len(matches)):
            # This is the starting position of the sequence
            start = matches[j] + 3

            # This is the ending position of the sequence
            end = matches[j] + len(seqs[i])

            bedS.write('chrR' + '\t' + str(start) + '\t' + str(end)  + '\t' + names[i] + '\n')

            # Very important, only generate a reference for the first fold (all other folds simulate CNV)
            if i == 0 and j == 0:
                bedR.write('chrR' + '\t' + str(start) + '\t' + str(end)  + '\t' + names[i] + '\n')

    sim = sim.replace('@', 'N')

    chrR.write('>chrR\n')
    chrR.write(ref + '\n')

    chrS.write('>chrR\n')
    chrS.write(sim + '\n')

    bedR.close()
    bedS.close()
    chrR.close()
    chrS.close()    

# This one works... We know...
def algo_2():
    
    ref = None
    
    # This is our reference genome!
    with open('CVA002.v015.fa', 'r') as f:
        ref = f.readlines()
    
    # We only care the sequence
    ref = ref[1]

    print ('Number of bases in the original: ' + str(len(ref)))
    
    # How many bases we want to take?
    n = 600
    
    # Where to start?
    i = 10673
    
    # How many copies we want?
    fold = 19
    
    j = i + n
    
    # This is our CNV
    copy = 'CGAATTGCCGCCAAGCAACGGTTCAAAGCATACTTAAGGCATCAGTTATTTTTATAGTGTAGTCGTTTAGATAAAACCGTCGCCTAATAGGTAAATGTGGGGAGGAAATTAAATGTCGAGGGAGGGTTGCGAATATGGATAATCCCTGCTTAACGACCTAAGAGTTTATTTCAGGGACTCAACGAAATTCCACTACTGAGATACCTGTTTTTAAAGCATCGACGTGGGATACAAAACCTTGTACTATCTTCTTCTCTTCAGGAGCCGTCTAACTGAATTAATAGTTAAATGCATCCTCTACATAGATGTGACCTGGATGCGCATTGGGTATCAGAAGTGGACCTATGGTATAAGTCTTTTGCACACGTTTTGGCTATCATAGGTTGTAGAATGCGACAGTGAGAAATTTTTTCCCCCCGTATGATCAAGAGCGGCGGTCTAACATCGTCATATCTAACATATGTTCATCTCCCATTACAGAAAAAATCTGGTGTGAAATAATGGATCTACACCTAACGTAGCCTCCGGAAGTACTGGGGCTTGAAAAGTTTTCACCAATACGAGCCAAATGAGAGGATGTTCAAGTCTTGGGAGTCCGCA'
    #copy = ref[i:j]
        
    # Starting for the copied sequence
    start = i
    
    # Ending fo the copied sequence
    end = j
    
    # Annotation for the reference
    bed1  = open('Algo2/chrT.bed', 'w')

    # Annotation for the sample
    bed2 = open('Algo2/chrF.bed', 'w')    

    # Where the copied sequence is (useful for visualizing alignments)
    bed1.write('chrT\t' + str(i) + '\t' + str(j) + '\t' + 'Copied\n')    

    # Assume we'll insert copies after the sequence...
    bed2.write('chrT\t' + str(i) + '\t' + str(j) + '\t' + 'Copied\n')        

    print ('Number of bases in CNV: ' + str(len(copy)))
    print ('The initial sequences in the copy is: ' + copy[:10])

    # Chromosome for the simulation
    chrF = open('Algo2/chrF.fa',  'w')

    for k in range(fold):
        ki = (2*j) + (k * 10000) # Separated by 10000...
        ref = ref[:ki] + '@' + copy + ref[ki:]

    i = 1
    for k in range(len(ref)):
        if ref[k] == '@':
            bed2.write('chrT\t' + str(k+1) + '\t' + str(k + n + 1) + '\t' + 'chrF_' + str(i) + '\n')
            i = i + 1

    ref = ref.replace('@', 'N')

    chrF.write('>chrF\n')
    chrF.write(ref + '\n')

    print ('Final number of bases in CNV: ' + str(len(ref)))

    bed1.close()
    bed2.close()
    chrF.close()
    chrT.close()


def generate(seqs, names, refs, chrC, bedA):

    # The position of the genome
    n = 0

    # How many useless paddings?
    paddings = 1000

    #
    # To simulate for CNV, we'll need to generate a genome with a copy of each sequin in the standard. But how?
    #

    for i in range(len(seqs)):
        
        if (names[i].endswith('_A')):
            fold = 1
        elif (names[i].endswith('_B')):
            fold = 2
        elif (names[i].endswith('_C')):
            fold = 4
        elif (names[i].endswith('_D')):
            fold = 8
        else:
            raise Exception('Unknown ' + names[i])
        
        for k in range(paddings):
            chrC.write('N')
        n = n + paddings
        
        # Repeat for all the copy times
        for j in range(1):
        #for j in range(fold):

            chrC.write(seqs[i])

            # Generate the corresponding annoation
            bedA.write('chrC\t' + str(n) + '\t' + str(n + len(seqs[i])) + '\t' + names[i] + '_' + str(j) + '\n')        
            n = n + len(seqs[i])

            for k in range(paddings):
                chrC.write('N')
            n = n + paddings

if __name__ == '__main__':
    
#    algo_2()
#    algo_3()
    #simulate()
    #simulate_3()
    #whereAreTheReads()
    #simulate_3()
    simulate_5()
    #generateGenomeFromUnjoined()
    #generateRandomChr21()

    #raise Exception('IT IS OK TO SEE ME!')

    #refs = []

    #with open('CVA002.v015.fa', 'r') as f:
    #    for line in f:
    #        refs.append(line) 

    #with open("unjoined.fa", "r") as f:
    #    names = []
    #    seqs  = []
    #    for line in f:
    #        if (line[0] == '>'):
    #            line = line.replace('>', '')                
    #            names.append(line.strip())
    #        else:
    #            seqs.append(line.strip())

    #chrC = open('chrC.fa', 'w')
    #bedA = open('chrC.bed', 'w')

    #chrC.write('>chrC\n')

    #for i in range(1):
    #    generate(seqs, names, refs, chrC, bedA)

    #chrC.write(refs[1])
    #chrC.write('\n')
    
    #bedA.close()
    #chrC.close()
