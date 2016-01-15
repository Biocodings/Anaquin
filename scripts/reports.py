#!/usr/bin/python

import os
import sys
import math

########################################################
#                                                      #
#         User options. Adjust them for yourself       #
#                                                      #
########################################################

ANAQUIN_PATH = '/Users/tedwong/Sources/QA/anaquin'


########################################################
#                                                      #
#        DON'T modify anything below this point!       #
#                                                      #
########################################################

EXPECT_NUM = 'Number'
EXPECT_STR = 'String'

# Eg: /home/tedwong/Sources/CountTable.txt
EXPECT_FILE = 'File'

# Eg: /home/tedwong/Sources/CountTable_1.txt,/home/tedwong/Sources/CountTable_2.txt
EXPECT_FILES = 'Files'

# Eg: ['1',2',3','4']
EXPECT_LIST = 'List'

# Where the temporary files are saved
TEMP_PATH = 'temp'

# Execute an Anaquin request
def anaquin(tool, args, config, needMixture=True):
    
    #
    # Generate a full Anaquin command. A full command would need the mixture and reference annotation.
    # However, not all tool would require, say, a mixture.
    #
    
    req = ANAQUIN_PATH + ' -t ' + tool + ' -m data/trans/MTR004.v013.csv -rgtf data/trans/ATR001.v032.gtf '
    
    # Now add up the arguments
    req = req + '-o ' + TEMP_PATH + ' ' + args

    print(req)
    #os.system('/Users/tedwong/Sources/QA/anaquin -o  ' + TEMP_PATH + ' ' + args)

###################################
#                                 #
#        Accessor functions       #
#                                 #
###################################

def checkFiles(root, files, max=None):

    # Turn the object into a list
    files = files.split(',')
    
    for i in range(0, len(files)):
        files[i] = root + files[i]
        if not(os.path.isfile(files[i])):
            raise Exception(files[i] + ': No such file')
            
    if (max is not None and len(files) > max):
        raise Exception('Too many files specified for ' + str(files))
            
    return ','.join(files)

def get(config, key, formats=None, optional=False):
    if (key in config):
        val = config[key]
        
        if (type(formats) is list and not(val in formats)):
            raise Exception('Key ' + key + ' must be ' + str(formats))
        elif (formats == EXPECT_FILE):
            return (checkFiles(root(config), val, 1))            
        elif (formats == EXPECT_FILES):
            return (checkFiles(root(config), val))
        return config[key]
    elif (optional == True):
        return None
    else:
        raise Exception('Failed to find ' + key + '. Please check and try again.')

# Return the root path
def root(config):
    path = get(config, 'ROOT_PATH', optional=True)
    
    # Nothing is given, assume everything is relative to the current directory
    if (path == None):
        return ''

    # Now, we can append a path afterward
    else:
        return path + '/'

def name(config):
    return get(config, 'REPORT_NAME', optional=True)
    
def mixture(config):
    return get(config, 'MIX_FILE', EXPECT_FILE)

# Return the column for mixture 1
def mixture_1(config):
    return get(config, 'MIX_1_COL', EXPECT_NUM)

# Return the column for mixture 2
def mixture_2(config):
    return get(config, 'MIX_2_COL', EXPECT_NUM, optional=True)


#################################
#                               #
#        Sequin functions       #
#                               #
#################################


# Create a report for TransQuin
def transQuin(config):
    
    metrics = get(config, 'DIFF_LEVEL', ['Gene', 'Isoform', 'Exon'])

    #############################################
    #                                           #
    #  1. Generating statistics for alignments  #
    #                                           #
    #############################################

    # Factors for the alignments
    factors = get(config, 'ALIGN_FACTORS')    

    # Alignment files
    alignFiles = get(config, 'ALIGN_FILE', EXPECT_FILES)

    #
    # Generate a request for TransQuin for differential analysis. For example:
    #
    #    anaquin TransAlign -factors 1,1,1,2,2,2 -ufiles C1.BAM,C2.BAM,C3.BAM
    #

    req = 'TransAlign -factors ' + factors + ' -ufiles ' + alignFiles
    
    # Execute the command
    anaquin('TransAlign', req, config)

    return

    ########################################################
    #                                                      #
    #  2. Generating statistics for differential analysis  #
    #                                                      #
    ########################################################

    # Factors for differential analysis
    factors = get(config, 'DIFF_FACTORS')
    
    countSoft  = get(config, 'COUNT_SOFT', { 'HTSeqCount' })    
    countFiles = get(config, 'DIFF_COUNT', EXPECT_FILES)

    print('Counting software: ' + countSoft)

    diffSoft = get(config, 'DIFF_SOFT', ['Cuffdiff', 'edgeR', 'DEXSeq2', 'DESeq2'])
    diffFile = get(config, 'DIFF_FILE', EXPECT_FILE)
    
    print('Differential software: ' + diffSoft)

    #
    # 5. Generate a request for TransQuin for differential analysis. For example:
    #
    #      anaquin TransDiff -soft DESEq2 -ufiles C1.txt,C2.txt,C3.txt
    #
    
    req = '-factors ' + factors + ' -countSoft ' + countSoft + ' -diffSoft ' + diffSoft + ' -countFiles ' + countFiles + ' -diffFile ' + diffFile

    # Execute the command
    anaquin('TransDiff', req, config)
    
def parse(file):

    print ('Parsing: ' + file)
    
    #
    # We simply need to parse all the key-value pairs. It's not the function's responsibility
    # to validate the inputs.
    #
    
    dict = {}

    for line in open(file):
        line = line.strip()
        
        if (len(line) == 0 or line[0] == '#'):
            continue
        
        toks = line.split('=')
        
        if (len(toks) != 2):
            raise Exception('Syntax error: ' + line)

        key = toks[0].strip()
        val = toks[1].strip()
        
        if (key in dict):
            dict[key] = dict[key] + ',' + val
        else:
            dict[key] = val

    print ('Parsing completed. ' + str(len(dict)) + ' keys found.')
    return dict

if __name__ == '__main__':

    mode = sys.argv[1]

    # The raw SAM file  
    file = sys.argv[2]

    if (mode == 'TransQuin'):
        transQuin(parse(file))
        
        
#        -t TransDiff -m data/trans/MTR004.v013.csv -rgtf data/trans/ATR001.v032.gtf -t TransDiff -factors 1,1,1,2,2,2 -soft DESEq2 -level gene -ufiles /Users/tedwong/Desktop/K_562/edgeR_Gene.csv 
        
        
#        /Users/tedwong/Sources/QA/anaquin   temp -factors -countSoft HTSeqCount -diffSoft edgeR -countFiles combined_counts_full/K_RMXA1v2.htseq.counts.combined,combined_counts_full/K_RMXA2v2.htseq.counts.combined,combined_counts_full/K_RMXA3v2.htseq.counts.combined,combined_counts_full/G_RMXB1v2.htseq.counts.combined,combined_counts_full/G_RMXB2v2.htseq.counts.combined,combined_counts_full/G_RMXB3v2.htseq.counts.combined -diffFile /Users/tedwong/Desktop/K_562/diffs.csv


#-t TransDiff -m data/trans/MTR004.v013.csv -rgtf data/trans/ATR001.v032.gtf -o temp TransAlign -factors 1,1,1,2,2,2 -ufiles /Users/tedwong/Desktop/K_562/Aligns/A1/accepted_hits.bam,/Users/tedwong/Desktop/K_562/Aligns/A2/accepted_hits.bam,/Users/tedwong/Desktop/K_562/Aligns/A3/accepted_hits.bam,/Users/tedwong/Desktop/K_562/Aligns/B1/accepted_hits.bam,/Users/tedwong/Desktop/K_562/Aligns/B2/accepted_hits.bam,/Users/tedwong/Desktop/K_562/Aligns/B3/accepted_hits.bam
