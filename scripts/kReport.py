#
# This file implements the k-mer quantification pipeline for Anaquin.
#

import io
import os
import sys
import tempfile
import subprocess

# Where Anaquin is located
ANAQUIN_PATH = 'anaquin'

# Where the temporary files are saved
TEMP_PATH = tempfile.gettempdir()

EXPECT_NUM = 'Number'
EXPECT_STR = 'String'

# Eg: /home/tedwong/Sources/CountTable.txt
EXPECT_FILE = 'File'

# Eg: /home/tedwong/Sources/CountTable_1.txt,/home/tedwong/Sources/CountTable_2.txt
EXPECT_FILES = 'Files'

# Eg: ['1',2',3','4']
EXPECT_LIST = 'List'



# Do we want to do unit testing?
__unitTesting__ = True

# Global variable
__mode__ = None

def download(url):
    path = tempfile.gettempdir()
    path = path + os.sep + str(uuid.uuid4())
    urllib.urlretrieve(url, path)    
    return path

def execute(cmd):
    print(cmd)
    os.system(cmd)

# Execute a FusQuin command
def rFusQuin(tool, args, config, onlyPrint=False):

    # Eg: /Users/tedwong/Desktop/K_562
    root = get(config, 'ROOT_PATH')

    mix   = appendRoot(root, get(config, 'MIX_FILE')) # CSV format
    cVars = appendRoot(root, get(config, 'CHRT_VAR')) # VCF format

    # Construct mixture and reference annoation
    req = ANAQUIN_PATH + os.sep + 'anaquin -t ' + tool + ' -m ' + mix + ' -rvcf ' + cVars

    # Now add up the arguments
    req = req + ' -o ' + TEMP_PATH + args

    print(req + '\n')
    
    # Execute the Anaquin request
    if not onlyPrint:
        os.system(req)

# Execute a VarQuin command
def rVarQuin(tool, args, config, onlyPrint=False):

    # Eg: /Users/tedwong/Desktop/K_562
    root = get(config, 'ROOT_PATH')

    mix   = appendRoot(root, get(config, 'MIX_FILE')) # CSV format
    cVars = appendRoot(root, get(config, 'CHRT_VAR')) # VCF format

    # Construct mixture and reference annoation
    req = ANAQUIN_PATH + os.sep + 'anaquin -t ' + tool + ' -m ' + mix + ' -rvcf ' + cVars

    # Now add up the arguments
    req = req + ' -o ' + TEMP_PATH + args

    print(req + '\n')
    
    # Execute the Anaquin request
    if not onlyPrint:
        os.system(req)

# Execute an TransQuin command
def rTransQuin(tool, args, config, onlyPrint=False, subdir=''):
    
    # Eg: /Users/tedwong/Desktop/K_562
    root = get(config, 'ROOT_PATH')

    names   = get(config, 'NAMES')
    factors = get(config, 'FACTORS')
    
    mix   = appendRoot(root, get(config, 'MIX_FILE'))   # CSV format
    c_ref = appendRoot(root, get(config, 'CHRT_ANNOT')) # GTF format
    e_ref = appendRoot(root, get(config, 'EXP_ANNOT'))  # GTF format
    
    names   = get(config, 'NAMES')    
    factors = get(config, 'FACTORS')
    
    if (len(subdir)):
        subdir = '/' + subdir
    
    # Construct mixture and reference annoation
    req = ANAQUIN_PATH + os.sep + 'anaquin -t ' + tool + ' -m ' + mix + ' -rgtf ' + c_ref + ' -rexp ' + e_ref
    
    # Now add up the arguments
    req = req + ' -o ' + TEMP_PATH + subdir + ' -factors ' + factors + ' -names ' + names + ' ' + args

    print(req + '\n')

    # Execute the Anaquin request
    if not onlyPrint:
        os.system(ANAQUIN_PATH + ' -o  ' + TEMP_PATH + ' ' + args)

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
    def writeRCode(file, output, src, title, nPlots, description, height=None):
        tmp = (tempfile.NamedTemporaryFile())
        if (output == 'RMarkdown'):
            file.write('\n## ' + title + '\n\n')
            file.write(description + '\n\n')
            file.write('```{r results=''\'hide\''', message=FALSE, warning=FALSE, echo=FALSE}\n')

            if (height == None):
                file.write('png(filename="' + tmp.name + '%01d")\n')
            else:
                file.write('png(filename="' + tmp.name + '%01d", height=' + str(height) + ' )\n')

            file.write('source("' + src + '")\n')
            file.write('dev.off()\n')
            file.write('```\n')

            for i in range(0, nPlots):
                file.write('\n![](' + tmp.name + str(i+1) + ')')

            file.write('\n')

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

    def addRCode(self, title, file, description, nPlots, height=None):
        self.items.append({ 'type': 'rCode', 'title': title, 'value': file, 'nPlots': nPlots, 'description': description, 'height': height })
        
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
                Language.writeRCode(file, output, TEMP_PATH + os.sep + item['value'], item['title'], item['nPlots'], item['description'], item['height'])
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
        
    def addRCode(self, title, file, description='', nPlots=1, height=None):
        self.current.addRCode(title, file, description, nPlots, height)
        
    def generate(self, file, output):
        print('\n-------------------------------------')
        print('Generating ' + file + ' for ' + output)
        
        if (output is 'RMarkdown'):
            file = open(file, 'w')

            header = """---
title: "Anaquin: %s"
header-includes: \usepackage{graphicx}
output: 
    pdf_document:
        keep_tex: true
        toc: yes
        toc_depth: 2
---""" % __mode__

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

def TransQuin(config, output):
    
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

    r.startChapter('Statistics (TransAlign)')

    for i in range(0, len(names)):
        r.addTextFile('Alignment statistics for: ' + names[i], names[i] + os.sep + 'TransAlign_summary.stats', )

    r.endChapter()

    ############################################
    #                                          #
    #  2. Generating statistics for assembly   #
    #                                          #
    ############################################

    print ('----------------------- Assembly -----------------------\n')

    # Assembly software
    soft = get(config, 'ASSEMBLY_SOFT', { 'Cufflinks', 'StringTie' })

    # Assembly files
    files = get(config, 'ASSEMBLY_FILE', EXPECT_FILES)
    
    #
    # Generate a request for TransAssembly for assembly analysis. For example:
    #
    #    anaquin TransAssembly -m ... -rgtf ... -factors 1,1,1,2,2,2 -names A1,A2,A3... -ufiles A1.gtf,A2.gtf,A3.gtf...
    #

    req = '-soft ' + soft + ' -ufiles ' + files
    
    # Execute the command
    anaquin('TransAssembly', req, config, onlyPrint=True)

    r.startChapter('Statistics (TransAssembly)')

    # Add summary statistics for each replicate
    for i in range(0, len(names)):
        r.addTextFile('Assembly statistics for: ' + names[i], names[i] + os.sep + 'TransAssembly_summary.stats', )

    r.endChapter()

    #############################################################
    #                                                           #
    #  3. Generating statistics for expression analysis (Gene)  #
    #                                                           #
    #############################################################

    print ('----------------------- Expression (Gene) -----------------------\n')

    soft  = get(config, 'EXP_G_SOFT', { 'Cufflinks', 'StringTie' })
    files = get(config, 'EXP_G_FILE', EXPECT_FILES)

    #
    # Generate a request for TransExp for expression analysis. For example:
    #
    #    anaquin TransExp -m ... -rgtf ... -lvl gene -factors 1,1,1,2,2,2 -names A1,A2,A3... -ufiles C1.exp,C2.exp,C3.exp...
    #

    req = '-soft ' + soft + ' -level gene  -ufiles ' + files
    
    # Execute the command for genes
    anaquin('TransExpress', req, config, onlyPrint=True, subdir='TG')

    r.startChapter('Statistics (Gene Expression)')

    dPooled  = 'The pooled scatter plot shows the expected abundance against measured abundance on the logarithm scale for all the replicates. This is done by plotting the average and standard deviation. The line is the fitted linear regression with 95% confidence interval drawn in black.\nA high R2 is desirable; the higher it is, the more accurate the gene expression experiment is.'
    dScatter = 'The scatter plot shows the relationship between expected abundance against measured abundance on the logarithm scale for a single replicate. The line is the fitted linear regression with 95% confidence interval drawn in black.\nA high R2 is desirable; the higher it is, the more accurate the gene expression experiment is.'

    r.addRCode('Pooled scatter plot for gene expression', 'TG' + os.sep + 'TransExpress_pooled.R', dPooled)
    r.addTextFile('Gene expression summary', 'TG' + os.sep + 'TransExpress_pooled.stats')

    # Add summary statistics for each replicate
    for i in range(0, len(names)):
        r.addTextFile('Gene expression statistics for: ' + names[i], 'TG' + os.sep + names[i] + os.sep + 'TransExpress_summary.stats', )
        r.addRCode('Gene expression scatter plot for: '  + names[i], 'TG' + os.sep + names[i] + os.sep + 'TransExpress_scatter.R', dScatter)

    r.endChapter()
    
    ################################################################
    #                                                              #
    #  4. Generating statistics for expression analysis (Isoform)  #
    #                                                              #
    ################################################################
    
    print ('----------------------- Expression (Isoform) -----------------------\n')

    # Expression software
    soft = get(config, 'EXP_I_SOFT', { 'Cufflinks', 'StringTie' })

    # Expression files
    files = get(config, 'EXP_I_FILE', EXPECT_FILES)

    dPooled  = dPooled.replace('gene',  'isoform')
    dScatter = dScatter.replace('gene', 'isoform')
    dMinor   = 'The Minor/Major plot shows the relative quantification of alternative spliced isoforms by measuring the minimum isoform as a fraction of the major isoform for each sequin gene. The accuracy of the quantification is typically dependent on sequence coverage and higher for high abundance genes. The concentration of the gene is shown by colors.'

    #
    # Generate a request for TransExp for expression analysis. For example:
    #
    #    anaquin TransExp -m ... -rgtf ... -lvl isoform -factors 1,1,1,2,2,2 -names A1,A2,A3... -ufiles C1.exp,C2.exp,C3.exp...
    #

    req = '-soft ' + soft + ' -level isoform  -ufiles ' + files
    
    # Execute the command for genes
    anaquin('TransExpress', req, config, onlyPrint=True, subdir='TI')

    r.startChapter('Statistics (Isoform Expression)')

    r.addRCode('Minor/Major plot',  'TI' + os.sep + 'TransExpress_Splice.R', dMinor)
    r.addRCode('Pooled scatter plot for isoform expression', 'TI' + os.sep + 'TransExpress_pooled.R', dPooled)
    r.addTextFile('Isoform expression summary', 'TI' + os.sep + 'TransExpress_pooled.stats')

    # Add summary statistics for each replicate
    for i in range(0, len(names)):
        r.addTextFile('Isoform expression statistics for: ' + names[i], 'TI' + os.sep + names[i] + os.sep + 'TransExpress_summary.stats')
        r.addRCode('Isoform expression scatter plot for: '  + names[i], 'TI' + os.sep + names[i] + os.sep + 'TransExpress_scatter.R', dScatter)

    r.endChapter()

    ########################################################
    #                                                      #
    #  5. Generating statistics for differential analysis  #
    #                                                      #
    ########################################################

    print ('----------------------- Differential -----------------------\n')

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

    r.startChapter('Statistics (Gene Expression Differential)')

    # Add summary statistics for each replicate
    r.addTextFile('Differential summary statistics', 'TransDiffs_summary.stats', )

    # Add scatter plot
    r.addRCode('Scatter plot', 'TransDiffs_scatter.R', 'The scatter plot shows the expected fold-change of genes against their measured fold-change. The more correlated they are the more accurate the experiment is. The R2 statistic measures how accurate the experiment is in detecting fold-changes. Accurate detection of fold-changes is crucial for differential analysis.\nThe accuracy of the quantification is typically dependent on abundance and higher for high abundance genes. The concentration of the gene is shown by colors.')

    # Add ROC plot
    r.addRCode('ROC plot', 'TransDiffs_ROC.R', 'When true differences in expression exist between samples in an experiment, those differences should be detected in differential expression tests; where no differences exist, no difference should be detected. The true-positive and true-negative sequins controls can be used in a receiver operator characteristic (ROC) curve analysis of rank ordered differential expression test P-values.')

    # Add MA plot
    r.addRCode('MA plot', 'TransDiffs_MA.R', 'The MA plot is used to study dependences between the log fold-ratio and the average normalized counts. The counts are normalized as they are adjusted for the library size. Typically, the variability is higher for less abundance genes.\n\nSequin data points are colored by log-ratio. They represent the mean log-fold changes for a given concentration level. Endogenous genes data points are shown as pink points.')

    # Add LODR plot
    r.addRCode('LODR plot', 'TransDiffs_LODR.R', 'Identifying differentially expressed genes is the objective of differential expression experiments; however, how much information is needed to have confidence that a given fold change in expression of transcripts will be detected? The LODR estimates can inform researchers of diagnostic power at a given fold change.\n\nAn LODR estimate for a particular fold change is the minimum signal, above which differentially expressed transcripts can be detected with a specified degree of confidence. The estimate for each ratio group is found based on the intersection of the model upper confidence interval upper bound with the P-value threshold.')

    r.endChapter()

    #########################################
    #                                       #
    #       6. Generating apprendix         #
    #                                       #
    #########################################

    r.startChapter('Apprendix: Sequin Alignment')

    for i in range(0, len(names)):
        r.addTextFile('Sequin statistics for: ' + names[i], names[i] + os.sep + 'TransAlign_quins.stats', )

    r.endChapter()

    # Generate a markup report (which can then be converted into various formats)
    r.generate('/Users/tedwong/Sources/QA/ABCD.RMarkdown', output)


def VarQuin(config, output):
    
    r = Report()
    
    # Eg: A1,A2,A3,B1,B2,B3
    names = getNames(config)

    #########################################
    #                                       #
    #       1. Generating alignments        #
    #                                       #
    #########################################
    
    print ('----------------------- Variant Alignment -----------------------\n')
    
    #
    # Generate a request for genome coverage. For example:
    #
    #    anaquin -t VarAlign -rbed data/VARQuin/AVA017.v032.bed -ufiles realigned.bam
    #

    files = get(config, 'ALIGN_FILE', EXPECT_FILES)

    req = ' -ufiles ' + files
    
    # Execute the command
    rVarQuin('VarAlign', req, config, onlyPrint=True)

    r.startChapter('Statistics (Genome Alignment)')

    for i in range(0, len(names)):
        r.addTextFile('Summary statistics for: ' + names[i], 'VarAlign_summary.stats', )

    r.endChapter()

    ########################################
    #                                      #
    #      2. Generating subsampling       #
    #                                      #
    ########################################

    print ('----------------------- Variant Subsampling -----------------------\n')

    #
    # Generate a request for genome coverage. For example:
    #
    #    anaquin -t VarSubsample -rbed data/VARQuin/AVA017.v032.bed -rexp chr21.bed -ufiles realigned.bam
    #

    files = get(config, 'COV_FILE', EXPECT_FILES)

    req = ' -ufiles ' + files
    
    # Execute the command
    rVarQuin('VarSubsampling', req, config, onlyPrint=True)

    r.startChapter('Statistics (Genome Subsampling)')

    for i in range(0, len(names)):
        r.addTextFile('Summary statistics for: ' + names[i], 'VarSubsample_summary.stats', )

    r.endChapter()

    #########################################
    #                                       #
    #    3. Generating variant discovery    #
    #                                       #
    #########################################

    print ('----------------------- Variant Discovery -----------------------\n')

    files = get(config, 'VAR_FILE', EXPECT_FILES)
    soft  = get(config, 'VAR_SOFT', { 'VarScan', 'GATK', 'FreeBayes' })

    # File name of the standards
    stand  = get(config, 'STAND_FILE', EXPECT_FILES)

    #
    # Generate a request for allele frequency. For example:
    #
    #    anaquin -t VarDiscover -rvcf data/VARQuin/AVA009.v032.vcf -rbed data/VARQuin/AVA017.v032.bed -m data/VARQuin/MVA011.v013.csv -soft VarScan -ufiles varscan.tab
    #
    
    req = ' -rbed ' + stand + ' -soft ' + soft + ' -ufiles ' + files
    
    # Execute the command
    rVarQuin('VarDiscover', req, config)

    r.startChapter('Statistics (Variant Discovery)')

    for i in range(0, len(names)):
        r.addTextFile('Summary statistics for: ' + names[i], 'VarDiscover_summary.stats', )
        r.addRCode('ROC Curve', 'VarDiscover_ROC.R', '', height=770)
        r.addRCode('LODR for SNP and Indel', 'VarDiscover_LODR.R', '')

    r.endChapter()

    ########################################
    #                                      #
    #    4. Generating allele frequency    #
    #                                      #
    ########################################

    print ('----------------------- Allele Frequency -----------------------\n')

    #
    # Generate a request for allele frequency. For example:
    #
    #    anaquin -t VarAllele -rvcf data/VARQuin/AVA009.v032.vcf -m data/VARQuin/MVA011.v013.csv -soft VarScan -ufiles varscan.tab 
    #

    req = ' -soft ' + soft + ' -ufiles ' + files
    
    # Execute the command
    rVarQuin('VarAllele', req, config)

    r.startChapter('Statistics (Allele Frequency)')

    for i in range(0, len(names)):
        r.addTextFile('Summary statistics for: ' + names[i], 'VarAllele_summary.stats', )
        r.addRCode('Scatter plot (SNPs + Indels)', 'VarAllele_scatter.R', '', nPlots=2) # TOOD: Fix this, should be 3 for indels
        r.addRCode('Coverage plot for SNP and Indel', 'VarAllele_coverage.R', '')

    r.endChapter()

    ##############################################
    #                                            #
    #           5. Generating Apprendix          #
    #                                            #
    ##############################################

    r.startChapter('Apprendix: Allele frequency')

    for i in range(0, len(names)):
        r.addTextFile('Allele frequency statistics for: ' + names[i], 'VarAllele_quins.csv', )

    r.endChapter()

    # Generate a markup report (which can then be converted into various formats)
    r.generate('report.Rmd', output)


def FusQuin(config, output):
    
    r = Report()
    
    # Eg: A1,A2,A3,B1,B2,B3
    names = getNames(config)

    ##############################################
    #                                            #
    #         1. Generating FusionDiscover       #
    #                                            #
    ##############################################

    print ('----------------------- Fusion Discovery -----------------------\n')

    #
    # Generate a request for fusion discovery. For example:
    #
    #    anaquin -t FusionDiscover -soft star -rfus2.ref -ufiles star-fusion.fusion_candidates.txt 
    #

    files = get(config, 'FUS_FILE', EXPECT_FILES)
    soft  = get(config, 'FUS_SOFT', { 'StarFusion', 'TopHatFusion', })
    req   = ' -soft ' + soft + ' -ufiles ' + files + ' -rfus data/fusion/AFU004.v032.ref '

    # Execute the command
    rFusQuin('FusionDiscover', req, config)

    r.startChapter('Statistics (Fusion Discovery)')

    for i in range(0, len(names)):
        r.addTextFile('Summary statistics for: ' + names[i], 'FusionDiscover_summary.stats', )

    r.endChapter()

    ##############################################
    #                                            #
    #         2. Generating FusionExpress        #
    #                                            #
    ##############################################

    print ('----------------------- Fusion Expression -----------------------\n')

    #
    # Generate a request for fusion expression. For example:
    #
    #    anaquin -t FusionExpress -soft star -rbed data/fusion/AFU006.v032.bed -m data/fusion/MFU007.v013.csv -ufiles SJ.out.tab
    #    anaquin -t FusionExpress -soft star -rbed data/fusion/AFU006.v032.bed -m data/fusion/MFU007.v013.csv -ufiles star-fusion.fusion_candidates.txt
    #

    req = ' -soft ' + soft + ' -ufiles ' + files + ' -rfus data/fusion/AFU004.v032.ref '

    # Execute the command
    rFusQuin('FusionExpress', req, config)

    r.startChapter('Statistics (Fusion Expression)')

    for i in range(0, len(names)):
        r.addTextFile('Summary statistics for: ' + names[i], 'FusionExpress_summary.stats', )

    r.endChapter()

    ##############################################
    #                                            #
    #   3. Generating for differential analysis  #
    #                                            #
    ##############################################

    print ('----------------------- Fusion Differential -----------------------\n')

    #
    # Generate a request for fusion differential. For example:
    #
    #    anaquin -t FusionDiff -soft star -rfus data/fusion/AFU004.v032.ref -rbed data/fusion/AFU006.v032.bed -m data/fusion/MFU007.v013.csv -utab SJ.out.tab -ufiles star-fusion.fusion_candidates.txt
    #

    req = ' -soft ' + soft + ' -rfus data/fusion/AFU004.v032.ref -rbed data/fusion/AFU006.v032.bed -m data/fusion/MFU007.v013.csv -utab SJ.out.tab -ufiles star-fusion.fusion_candidates.txt '

    # Execute the command
    rFusQuin('FusionDiff', req, config)

    r.startChapter('Statistics (Fusion Differential)')

    for i in range(0, len(names)):
        r.addTextFile('Summary statistics for: ' + names[i], 'FusionDiff_summary.stats', )

    r.endChapter()

    ##############################################
    #                                            #
    #           4. Generating Apprendix          #
    #                                            #
    ##############################################

    r.startChapter('Apprendix: Fusion Discover')

    for i in range(0, len(names)):
        r.addTextFile('Fusion Discover statistics for: ' + names[i], 'FusionDiscover_labels.csv', )

    r.endChapter()

    # Generate a markup report (which can then be converted into various formats)
    r.generate('report.RMarkdown', output)


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

def report2PDF(r, name, path):
    
    global TEMP_PATH
    TEMP_PATH = path
    
    # Generate markup documentation
    report = path + '/report.Rmd'

    # The script required to construct a PDF from the markdown
    r2pdf = path + '/r2pdf.R'
    
    r.generate(report, 'RMarkdown')

    f = open(r2pdf, 'w')
    f.write('library(Anaquin)\n')
    f.write('library(rmarkdown)\n')
    f.write('render("report.Rmd", "pdf_document")\n')
    f.close()
    
    # The script should be run relative to the output directory
    os.chdir(path)

    # Create the PDF report
    execute('R CMD BATCH ' + r2pdf)

    # Move it to where it's supposed to be
    execute('mv report.pdf ' + name + ' 2>/dev/null')

    #execute('rm ' + r2pdf)
    #execute('rm ' + report)

    print('PDF generated. Please check: ' + name)

# Running threads
thds = []

def addCmd(target, args):

    global thds    

    thd = threading.Thread(target=target, args=args)
    thds.append(thd)

def wait():
    for thd in thds:
        thd.start()

    for thd in thds:
        thd.join()

# Generate a VarQuin report with alignment-free k-mers
def VarQuinKM(anaq, path, index, mix, file1, file2):
    
    global TEMP_PATH
    TEMP_PATH = path
    
    r = Report()
    
    def VarKExpress(r):
        
        #########################################
        #                                       #
        #       1. Generating VarKExpress       #
        #                                       #
        #########################################
    
        # Eg: anaquin -t VarKExpress -m MVA011.v013.csv -rind AVA010.v032.index -ufiles LVA086.1_val_1.fq -ufiles LVA086.2_val_2.fq
        execute(anaq + ' -o ' + path + ' -t VarKExpress -soft kallisto -m ' + mix + ' -rind ' + index + ' -ufiles ' + file1 + ' -ufiles ' + file2)
    
        r.startChapter('Statistics (Expression)')
    
        r.addTextFile('Summary statistics', 'VarKExpress_summary.stats', )
        r.addRCode('Expected abundance vs measured abundance', 'VarKExpress_express.R', '')

        r.endChapter()
    
    def VarKAllele(r):
            
        #########################################
        #                                       #
        #       2. Generating VarKAllele        #
        #                                       #
        #########################################
    
        # Eg: anaquin -t VarKAllele -m MVA011.v013.csv -rind AVA010.v032.index -ufiles LVA086.1_val_1.fq -ufiles LVA086.2_val_2.fq 
        execute(anaq + ' -o ' + path + ' -t VarKAllele -soft kallisto -m ' + mix + ' -rind ' + index + ' -ufiles ' + file1 + ' -ufiles ' + file2)
    
        r.startChapter('Statistics (Allele Frequency)')
    
        r.addTextFile('Summary statistics', 'VarKAllele_summary.stats', )
        r.addRCode('Expected allele frequency vs measured allele frequency', 'VarKAllele_allele.R', '')

        r.endChapter()

    addCmd(VarKAllele,  [r])    
    addCmd(VarKExpress, [r])
    wait()

    ##############################################
    #                                            #
    #           3. Generating Apprendix          #
    #                                            #
    ##############################################

    r.startChapter('Apprendix:')
    r.addTextFile('Statistics for expression: ', 'VarKExpress_quins.csv', )
    r.endChapter()

    r.startChapter('Apprendix:')
    r.addTextFile('Statistics for allele frequency: ', 'VarKAllele_quins.csv', )
    r.endChapter()

    generatePDF(r, path, 'VarReport_report.pdf')

# Generate a TransQuin report with alignment-free k-mers
def TransQuinKM(anaq, path, index, mix, file1, file2):

    global TEMP_PATH
    TEMP_PATH = path

    r = Report()
    
    # In TransQuin, the file is the metadata
    file = file1

    def TransKExpress(r):
        
        ###########################################
        #                                         #
        #       1. Generating TransKExpress       #
        #                                         #
        ###########################################
        
        # Eg: anaquin -t TransKExpress -soft kallisto -rind data/TransQuin/ATR003.v032.index -m data/TransQuin/MTR002.v013.csv  -ufiles kallisto/abundance.tsv
        execute(anaq + ' -o ' + path + ' -t TransKExpress -soft kallisto -m ' + mix + ' -rind ' + index + ' -ufiles ' + file1 + ' -ufiles ' + file2)
    
        r.startChapter('Statistics (Expression)')
    
        r.addTextFile('Summary statistics', 'TransKExpress_summary.stats', )
        r.addRCode('Expected abundance vs measured abundance', 'TransKExpress_express.R', '')

        r.endChapter()

    def TransKDiff(r, meta):
    
        #########################################
        #                                       #
        #       2. Generating TransKDiff        #
        #                                       #
        #########################################
    
        # Eg: anaquin -t TransKDiff -rind data/TransQuin/ATR003.v032.index -m data/TransQuin/MTR004.v013.csv -soft sleuth -ufiles meta.csv
        execute(anaq + ' -o ' + path + ' -t TransKDiff -soft kallisto -m ' + mix + ' -rind ' + index + ' -ufiles ' + meta)
    
        r.startChapter('Statistics (Differential analysis)')
    
        r.addTextFile('Summary statistics', 'TransKDiff_summary.stats', )
        r.addRCode('Expected log-fold vs measured log-fold', 'TransKDiff_fold.R', '')
        #r.addRCode('ROC Curve', 'TransKDiff_ROC.R', '')

        r.endChapter()

    addCmd(TransKDiff, [r, file])
    wait()
    
    ##############################################
    #                                            #
    #           3. Generating Apprendix          #
    #                                            #
    ##############################################

    r.startChapter('Apprendix:')
    r.addTextFile('Statistics for differential: ', 'TransKDiff_quins.csv', )
    r.endChapter()

    generatePDF(r, path, 'TransReport_report.pdf')

# Create a PDF report for RnaQuin
def createRNReport(name, inputs):    
    r = Report()
    r.startChapter('RNAQuin Report')

    #
    # Expression analysis at isoform and gene level
    # Differential analysis at isoform and gene level
    #

    expIF = ['ExpressI/RnaExpress_summary.stats', 'ExpressI/RnaExpress_sequins.csv', 'ExpressI/RnaExpress_linear.R']
    expID = ['Summary Statistics (Isoform)', 'Detailed Statistics', 'Gene Expression']    
    expGF = ['ExpressG/RnaExpress_summary.stats', 'ExpressG/RnaExpress_sequins.csv', 'ExpressG/RnaExpress_linear.R']
    expGD = ['Summary Statistics (Gene)', 'Detailed Statistics', 'Gene Expression']

    diffIF = ['FoldI/RnaFold_summary.stats', 'FoldI/RnaFold_sequins.csv', 'FoldI/RnaFoldChange_fold.R']
    diffID = ['Summary Statistics (Isoform)', 'Detailed Statistics', 'Fold Plot (Isoform)']    
    diffGF = ['FoldG/RnaFold_summary.stats', 'FoldG/RnaFold_sequins.csv', 'FoldG/RnaFoldChange_fold.R']
    diffGD = ['Summary Statistics (Gene)', 'Detailed Statistics', 'Fold Plot (Gene)']
    
    files = expIF + expGF + diffIF + diffID
    descs = expID + expGD + diffGF + diffGD

    for i in range(0, len(files)):
        file = files[i]
        desc = descs[i]
        
        if (desc == ''):
            desc = file
        
        if file.endswith('R'):
            r.addRCode(desc, file, '')
        elif file.endswith('stats') and ('quins' in file):
            r.addTextFile(desc, file)
        else:
            r.addTextFile(desc, file)
            
    r.endChapter()
   
    report2PDF(r, name, inputs)

def run(cmd):
    try:
        print cmd
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        return e.output

def createText(file, txt):
    with io.FileIO(file, "w") as file:
        file.write(txt)
        
def runRScript(file):
    return run('Rscript ' + file)

# Check if an executable is installed (eg: Kallisto)
def checkInstall(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
    
def checkRInstall(pack):
    tmp = tempfile.gettempdir() + '/RPackInstalled.R'    
    cmd = pack + ' %in% rownames(installed.packages())'

    # Create an R-script for checking the existance of the package
    run('echo "require(' + pack + ')" > ' + tmp)

    # Eg: object 'ggplot2s' not found
    out = run('Rscript ' + tmp)
    
    return not ('no package' in out)

def RnaQuin(index, output, A1, A2, B1, B2):
    if checkInstall('kallisto') is None:
        raise Exception('Kallisto is not installed. Please consult the user guide on www.sequin.xyz and try again.')

    if checkInstall('R') is None:
        raise Exception('R is not installed. Please consult the user guide on www.sequin.xyz and try again.')

    if not checkRInstall('sleuth'):
        raise Exception('Sleuth R-package is not installed. Please consult the user guide on www.sequin.xyz and try again.')

    samps = []
    conds = []
    
    # For each replicate in the first sample...
    for i in range(len(A1)):
        samps.append('A' + str(i+1))
        conds.append('0')

        # Command for Kallisto
        quant = 'kallisto quant -b 500 -i ' + index + ' -o ' + (output + '/A' + str(i+1)) + ' ' + A1[i] + ' ' + A2[i]

        # Quantify the replicate
        run(quant)

    # For each replicate in the second sample...
    for i in range(len(B1)):
        samps.append('B' + str(i+1))
        conds.append('1')        

        # Command for Kallisto
        quant = 'kallisto quant -b 500 -i ' + index + ' -o ' + (output + '/B' + str(i+1)) + ' ' + B1[i] + ' ' + B2[i]

        # Quantify the replicate
        run(quant)
    
    print 'Kallisto quantification completed'
    
    #
    # Reference: https://www.biostars.org/p/143458/#148465
    #            https://www.biostars.org/p/157240/#157242
    #

    # Do we have at least two replicates in each sample?
    if len(A1) > 0 and len(B1) > 0:
        sleuth = """
                 library(sleuth)
    
                 # Where the Kallisto files are
                 path <- '%s'
    
                 # List of sample IDs
                 sample_id <- dir(file.path(path))
                 
                 samps <- c(%s)
                 conds <- c(%s)
    
                 # Construct full path for the samples
                 path <- paste(path, samps, sep='/')
    
                 s2c <- data.frame(sample=samps, condition=conds)
                 s2c <- dplyr::mutate(s2c, path=path)
    
                 so <- sleuth_prep(s2c, ~condition)
                 so <- sleuth_fit(so)
                 so <- sleuth_wt(so, 'condition1')
    
                 results <- sleuth_results(so, 'condition1')
    
                 write.csv(results, file='%s', row.names=FALSE, quote=FALSE)"""
                 
        samps = ["'" + s + "'" for s in samps]
        conds = ["'" + s + "'" for s in conds]        
    
        sleuth = sleuth % (output, ','.join(samps), ','.join(conds), output + '/sleuth.csv')
    
        # Create a R script for Sleuth
        createText(output + '/sleuth.R', sleuth)
        
        print output + '/sleuth.R'
    
        # Run Sleuth
        runRScript(output + '/sleuth.R')
        
        print 'Sleuth quantification completed'

#
# VarQuin currently supports Kallisto and Salmon.    
#
def VarKReport(file):
    pass

def parseSamples(argv, A1, A2, B1, B2):
    #
    # It's common to sample replicates for Rna-Seq. Parse those samples.
    #
        
    # What's the design for A? Eg: 3,3
    A = 2*int(sys.argv[4].split(',')[0]) # Paired-end
    
    # What's the design for B? Eg: 3,3
    B = 2*int(sys.argv[4].split(',')[1]) # Paired-end

    for i in sys.argv[5:5+A]:
        if len(A1) == len(A2):
            A1.append(i)
        else:
            A2.append(i)

    if B > 0:
        for i in sys.argv[5+A:5+A+B]:
            if len(B1) == len(B2):
                B1.append(i)
            else:
                B2.append(i)

#
# Eg: python ~/Sources/QA/scripts/kReport.py RnaQuin ARN024_v001.index /tmp/kallisto 3,3 LRN087.1_val_1.fq LRN087.2_val_2.fq LRN088.1_val_1.fq LRN088.2_val_2.fq LRN089.1_val_1.fq LRN089.2_val_2.fq LRN090.1_val_1.fq LRN090.2_val_2.fq LRN091.1_val_1.fq LRN091.2_val_2.fq LRN092.1_val_1.fq LRN092.2_val_2.fq
#
# Eg: python kReport.py RnaReport RnaKReport_report.pdf /tmp/kallisto
#
# Note: please give full path
#
if __name__ == '__main__':

    # Which sequin?
    __mode__ = sys.argv[1]

    # Eg: A.R.3.index
    file = sys.argv[2]

    # Where the outputs should be written to
    output = sys.argv[3]
    
    if not os.path.exists(output):
        os.makedirs(output)    

    if __mode__ == 'RnaQuin':
        
        A1 = []
        A2 = []
        B1 = []
        B2 = []

        parseSamples(sys.argv, A1, A2, B1, B2)
        RnaQuin(file, output, A1, A2, B1, B2)

    elif __mode__ == 'VarQuin':
        VarQuin(file, output, pair1, pair2)

    elif __mode__ == 'RnaReport':        
        createRNReport(file, output)

    elif __mode__ == 'VarReport':
        pass