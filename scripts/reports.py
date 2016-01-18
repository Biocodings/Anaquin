#!/usr/bin/python

import os
import sys
import math

########################################################
#                                                      #
#         User options. Adjust them for yourself       #
#                                                      #
########################################################

# Where Anaquin is located
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
    
    names   = get(config, 'NAMES')
    factors = get(config, 'FACTORS')
    
    #
    # Generate a full Anaquin command. A full command would need the mixture and reference annotation.
    # However, not all tool would require, say, a mixture.
    #
    
    # Construct mixture and reference annoation
    req = ANAQUIN_PATH + ' -t ' + tool + ' -m data/trans/MTR004.v013.csv -rgtf data/trans/ATR001.v032.gtf '
    
    # Now add up the arguments
    req = req + '-o ' + TEMP_PATH + ' ' + args

    print(req)
    
    # Execute the Anaquin request
    #os.system(ANAQUIN_PATH + ' -o  ' + TEMP_PATH + ' ' + args)

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
#        Report functions       #
#                               #
#################################

class Report:
    def __init__(self, output):
        self.output = output
        
    def addSummary(self):
        pass
        
    def addImage(self):
        pass
        
    def startSection():
        pass
        
    def endSection():
        pass

def createReport(file, report):
    pass

#################################
#                               #
#        Sequin functions       #
#                               #
#################################


# Create a report for TransQuin
def transQuin(config, factors, names, output):
    
    ###########################################
    #                                         #
    #    1. Table of contents                 #
    #    2. Alignments                        #
    #    3. Assembly                          #
    #    4. Gene expression                   #
    #    5. Isoform expression                #
    #    6. Alternative splicing expression   #
    #    7. Differential expression           #
    #                                         #
    ###########################################

    r = Report(output)
    
    #############################################
    #                                           #
    #  1. Generating statistics for alignments  #
    #                                           #
    #############################################

    # Alignment files
    alignFiles = get(config, 'ALIGN_FILE', EXPECT_FILES)

    #
    # Generate a request for TransQuin for differential analysis. For example:
    #
    #    anaquin TransAlign -m ... -rgtf ... -factors 1,1,1,2,2,2 -ufiles C1.BAM,C2.BAM,C3.BAM
    #

    req = 'TransAlign -factors ' + factors + ' -ufiles ' + alignFiles
    
    # Execute the command
    anaquin('TransAlign', req, config)


    ######################################################
    #                                                    #
    #  2. Generating statistics for expression analysis  #
    #                                                    #
    ######################################################

    # Expression files
    expSoft = get(config, 'EXP_SOFT', { 'Cuffdiff', 'StringTie' })

    #
    # Generate a request for TransExp for expression analysis. For example:
    #
    #    anaquin TransExp -m ... -rgtf ... -factors 1,1,1,2,2,2 -names A1,A2,A3... -ufiles C1.exp,C2.exp,C3.exp...
    #

    req = '-ufiles ' + expSoft
    
    # Execute the command
    anaquin('TransExp', req, config)
    
    
    #####################################################
    #                                                   #
    #  Generating statistics for differential analysis  #
    #                                                   #
    #####################################################

    metrics = get(config, 'DIFF_LEVEL', ['Gene', 'Isoform', 'Exon'])

    countSoft  = get(config, 'COUNT_SOFT', { 'HTSeqCount' })    
    countFiles = get(config, 'DIFF_COUNT', EXPECT_FILES)

    print('Counting software: ' + countSoft)

    diffSoft = get(config, 'DIFF_SOFT', ['Cuffdiff', 'edgeR', 'DEXSeq2', 'DESeq2'])
    diffFile = get(config, 'DIFF_FILE', EXPECT_FILE)
    
    print('Differential software: ' + diffSoft)

    #
    # Generate a request for TransQuin for differential analysis. For example:
    #
    #   anaquin TransDiff -m ... -rgtf ... -soft DESEq2 -ufiles P1.txt,P2.txt,P3.txt... -cfiles C1.txt,C2.txt,C3.txt...
    #
    
    req = '-factors ' + factors + ' -countSoft ' + countSoft + ' -diffSoft ' + diffSoft + ' -countFiles ' + countFiles + ' -diffFile ' + diffFile

    # Execute the command
    anaquin('TransDiff', req, config)


#################################
#                               #
#         Misc functions        #
#                               #
#################################
    
    
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
    
    output = 'RMarkdown'

    if (mode == 'TransQuin'):
        transQuin(parse(file), output)
        
