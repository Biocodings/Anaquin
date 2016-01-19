#!/usr/bin/python

#
# This script implements TransReport, VarReport, FusionReport and MetaReport.
#

import os
import sys
import math
import uuid


########################################################
#                                                      #
#         User options. Adjust them for yourself       #
#                                                      #
########################################################

# Where Anaquin is located
ANAQUIN_PATH = '/Users/tedwong/Sources/QA'


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
TEMP_PATH = ANAQUIN_PATH + os.sep + '__temp__' #str(uuid.uuid4())


# Do we want to do unit testing?
__unitTesting__ = True


# Execute an Anaquin request
def anaquin(tool, args, config, needMixture=True, onlyPrint=False):
    
    # Eg: /Users/tedwong/Desktop/K_562
    root = get(config, 'ROOT_PATH')

    names   = get(config, 'NAMES')
    factors = get(config, 'FACTORS')
    
    mix = appendRoot(root, get(config, 'MIX_FILE'))
    ref = appendRoot(root, get(config, 'REF_ANNOT'))
    
    names   = get(config, 'NAMES')    
    factors = get(config, 'FACTORS')
    
    #
    # Generate a full Anaquin command. A full command would need the mixture and reference annotation.
    # However, not all tool would require, say, a mixture.
    #
    
    # Construct mixture and reference annoation
    req = ANAQUIN_PATH + os.sep + 'anaquin -t ' + tool + ' -m ' + mix + ' -rgtf ' + ref
    
    # Now add up the arguments
    req = req + ' -o ' + TEMP_PATH + ' -factors ' + factors + ' -names ' + names + ' ' + args

    print(req + '\n')

    # Execute the Anaquin request
    #if not __unitTesting__ or onlyPrint:
    #    os.system(ANAQUIN_PATH + ' -o  ' + TEMP_PATH + ' ' + args)

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

def appendRoot(root, path):
    if root is None:
        return path
    else:
        return root + os.sep + path

def reportName(config):
    return get(config, 'REPORT_NAME', optional=True)
    
def mixture(config):
    return get(config, 'MIX_FILE', EXPECT_FILE)

# Eg: A1,A2,A3,B1,B2,B3
def getNames(config):
    return get(config, 'NAMES').split(',')


#################################
#                               #
#        Report functions       #
#                               #
#################################

#
# Supported markup languages:
#
#    - RMarkdown
#
class Language:
    @staticmethod    
    def writePage(file, output):
        if (output == 'RMarkdown'):
            file.write('\n\pagebreak\n\n')

    @staticmethod
    def writeText(file, output, text):
        if (output == 'RMarkdown'):
            file.write(text)

    @staticmethod
    def writeTextFile(file, output, src, title):
        with open(src, 'r') as src:
            text = src.read()        
        
        if (output == 'RMarkdown'):
            file.write('\n## ' + title + '\n\n')
            file.write('```{ eval=FALSE}\n')
            file.write(text)
            file.write('\n```\n\n')

    @staticmethod
    def writeRCode(file, output, src, title):
        with open(src, 'r') as src:
            text = src.read()        
        
        if (output == 'RMarkdown'):
            file.write('\n## ' + title + '\n\n')
            file.write('```{r results=''\'hide\''', message=FALSE, warning=FALSE, echo=FALSE}\n')
            file.write(text)
            file.write('\n```\n\n')

class Chapter:
    def __init__(self, title):
        self.items = []
        self.title = title

    def addPage(self):
        self.items.append({ 'type': 'page', 'value': None })
    
    def addTextFile(self, title, file):
        self.items.append({ 'type': 'textFile', 'title': title, 'value': file })
        
    def addImage(self, title, file):
        self.items.append({ 'type': 'image', 'title': title, 'value': file })

    def addRCode(self, title, file):
        self.items.append({ 'type': 'rCode', 'title': title, 'value': file })
        
    def generate(self, file, output):
        Language.writePage(file, output)
        Language.writeText(file, output, '\n# ' + self.title + '\n')
        #Language.writeText(file, output, '---\n\n')

        for i in range(0, len(self.items)):
            item = self.items[i]
            
            if item['type'] == 'page':
                Language.writePage(file, output)                
            elif item['type'] == 'textFile':
                Language.writeTextFile(file, output, TEMP_PATH + os.sep + item['value'], item['title'])
                Language.writePage(file, output)
            elif item['type'] == 'rCode':
                Language.writeRCode(file, output, TEMP_PATH + os.sep + item['value'], item['title'])
                Language.writePage(file, output)                
            elif item['type'] == 'image':
                Language.writeImage(file, output, TEMP_PATH + os.sep + item['value'], item['title'])
                Language.writePage(file, output)                
            else:
                raise Exception('Unknown item: ' + str(item))
        
class Report:
    def __init__(self):
        self.chapters = []
    
    def startChapter(self, title):
        self.current = Chapter(title)
        
    def endChapter(self):
        self.chapters.append(self.current)
        self.current = None
        
    def addPage(self):
        self.current.addPage()

    def addTextFile(self, title, file):
        self.current.addTextFile(title, file)
        
    def addImage(self, title, file):
        self.current.addImage(title, file)
        
    def addRCode(self, title, file):
        self.current.addRCode(title, file)
        
    def generate(self, file, output):
        print('\n-------------------------------------')
        print('Generating ' + file + ' for ' + output)
        
        if (output is 'RMarkdown'):
            file = open(file, 'w')

            header = """---
title: "Anaquin: TransQuin Report"
header-includes: \usepackage{graphicx}
output: 
    pdf_document:
        keep_tex: true
        toc: yes
        toc_depth: 2
---"""

            file.write(header + '\n\n')

            for i in range(0, len(self.chapters)):
                self.chapters[i].generate(file, output)

def createReport(file, report):
    pass


#################################
#                               #
#        Sequin functions       #
#                               #
#################################


# Create a report for TransQuin
def transQuin(config, output):
    
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

    r = Report()
    
    # Eg: A1,A2,A3,B1,B2,B3
    names = getNames(config)
    
    #############################################
    #                                           #
    #  1. Generating statistics for alignments  #
    #                                           #
    #############################################

    print ('----------------------- Alignments -----------------------\n')

    # Alignment files
    files = get(config, 'ALIGN_FILE', EXPECT_FILES)

    #
    # Generate a request for TransQuin for differential analysis. For example:
    #
    #    anaquin TransAlign -m ... -rgtf ... -factors 1,1,1... -names A1,A2,A3... -ufiles C1.BAM,C2.BAM,C3.BAM
    #

    req = '-ufiles ' + files
    
    # Execute the command
    anaquin('TransAlign', req, config, onlyPrint=True)

    r.startChapter('TransQuin Alignment')

    # Add summary statistics for each replicate
    for i in range(0, len(names)):
        r.addTextFile('Alignment summary statistics for: ' + names[i], names[i] + os.sep + 'TransAlign_summary.stats', )
        
    # Add sequin statistics for each replicate
    for i in range(0, len(names)):
        r.addTextFile('Alignment sequin statistics for: ' + names[i], names[i] + os.sep + 'TransAlign_quins.stats', )

    r.endChapter()


    ############################################
    #                                          #
    #  2. Generating statistics for assembly   #
    #                                          #
    ############################################

    print ('----------------------- Assembly -----------------------\n')

    # Expression software
    soft = get(config, 'ASSEMBLY_SOFT', { 'Cufflinks', 'StringTie' })

    # Expression files
    #files = get(config, 'ASSEMBLY_FILE', EXPECT_FILES)
    
    
    ######################################################
    #                                                    #
    #  3. Generating statistics for expression analysis  #
    #                                                    #
    ######################################################

    print ('----------------------- Expression -----------------------\n')

    # Expression software
    soft = get(config, 'EXP_SOFT', { 'Cufflinks', 'StringTie' })

    # Expression files
    files = get(config, 'EXP_FILE', EXPECT_FILES)
    
    #
    # Generate a request for TransExp for expression analysis. For example:
    #
    #    anaquin TransExp -m ... -rgtf ... -factors 1,1,1,2,2,2 -names A1,A2,A3... -ufiles C1.exp,C2.exp,C3.exp...
    #

    req = '-soft ' + soft + ' -ufiles ' + files
    
    # Execute the command
    anaquin('TransExpress', req, config, onlyPrint=True)
    
    r.startChapter('TransQuin Expression')

    # Add summary statistics for each replicate
    for i in range(0, len(names)):
        r.addTextFile('Expression summary statistics for: ' + names[i], names[i] + os.sep + 'TransExpress_summary.stats', )
        
    # Add scatter plot
    for i in range(0, len(names)):
        r.addRCode('Expression scatter plot for: ' + names[i], names[i] + os.sep + 'TransExpress_scatter.R', )

    r.endChapter()

    
    #####################################################
    #                                                   #
    #  Generating statistics for differential analysis  #
    #                                                   #
    #####################################################

    lvl = get(config, 'DIFF_LEVEL', ['Gene', 'Isoform', 'Exon'])

    cSoft  = get(config, 'COUNT_SOFT', { 'HTSeqCount' })    
    cFiles = get(config, 'DIFF_COUNT', EXPECT_FILES)

    print('Counting software: ' + cSoft)

    dSoft = get(config, 'DIFF_SOFT', ['Cuffdiff', 'edgeR', 'DEXSeq2', 'DESeq2'])
    dFile = get(config, 'DIFF_FILE', EXPECT_FILE)
    
    print('Differential software: ' + dSoft)

    #
    # Generate a request for TransQuin for differential analysis. For example:
    #
    #   anaquin TransDiff -m ... -rgtf ... -soft DESEq2 -ufiles P1.txt,P2.txt,P3.txt... -cfiles C1.txt,C2.txt,C3.txt...
    #
    
    req = ' -soft ' + dSoft + ' -csoft ' + cSoft + ' -level ' + lvl + ' -cfiles ' + cFiles + ' -ufiles ' + dFile

    # Execute the command
    anaquin('TransDiff', req, config, onlyPrint=True)

    r.startChapter('TransQuin Differential')

    # Add summary statistics for each replicate
    r.addTextFile('Differential summary statistics', 'TransDiffs_summary.stats', )

    # Add ROC plots
    r.addRCode('ROC plot', 'TransDiffs_ROC.R', )

    # Add ROC plots
    r.addRCode('MA plot', 'TransDiffs_MA.R', )

    r.endChapter()

    # Generate a markup report (which can then be converted into various formats)
    r.generate('/Users/tedwong/Sources/QA/ABCD.RMarkdown', output)
    

#################################
#                               #
#       Parsing functions       #
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
        
