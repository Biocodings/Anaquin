#
# This file implements the k-mer quantification pipeline for Anaquin.
#

import io
import os
import sys
import numpy
import pandas
import tempfile
import subprocess
from multiprocessing.dummy import Pool as ThreadPool 

# Where the temporary files are saved
TEMP_PATH = tempfile.gettempdir()

# Eg: /home/tedwong/Sources/CountTable.txt
EXPECT_FILE = 'File'

# Eg: /home/tedwong/Sources/CountTable_1.txt,/home/tedwong/Sources/CountTable_2.txt
EXPECT_FILES = 'Files'

# Global variable
__mode__ = None

def execute(cmd):
    print(cmd)
    os.system(cmd)

def checkRInstall(pack):
    tmp = tempfile.gettempdir() + '/RPackInstalled.R'    
    cmd = pack + ' %in% rownames(installed.packages())'

    # Create an R-script for checking the existance of the package
    run('echo "require(' + pack + ')" > ' + tmp)

    # Eg: object 'ggplot2s' not found
    out = run('Rscript ' + tmp)
    
    return not ('no package' in out)

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

    diffIF = ['FoldI/RnaFoldChange_summary.stats', 'FoldI/RnaFoldChange_sequins.csv', 'FoldI/RnaFoldChange_fold.R']
    diffID = ['Summary Statistics (Isoform)', 'Detailed Statistics', 'Fold Plot (Isoform)']    
    diffGF = ['FoldG/RnaFoldChange_summary.stats', 'FoldG/RnaFoldChange_sequins.csv', 'FoldG/RnaFoldChange_fold.R']
    diffGD = ['Summary Statistics (Gene)', 'Detailed Statistics', 'Fold Plot (Gene)']
    
    files = expIF + expGF + diffIF + diffGF
    descs = expID + expGD + diffID + diffGD

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
    
def RnaQuin(index, output, exp):
    if checkInstall('kallisto') is None:
        raise Exception('Kallisto is not installed. Please consult the user guide on www.sequin.xyz and try again.')
    elif checkInstall('R') is None:
        raise Exception('R is not installed. Please consult the user guide on www.sequin.xyz and try again.')
    elif not checkRInstall('sleuth'):
        raise Exception('Sleuth R-package is not installed. Please consult the user guide on www.sequin.xyz and try again.')

    samps = []
    conds = []
    
    def runKallisto(cmd):
        # Command for Kallisto
        quant = 'kallisto quant -b 500 -i ' + index + ' -o ' + cmd

        # Quantify the replicate
        run(quant)

    cmds = []

    # For each replicate for mixture A...
    for i in range(len(exp.A)):
        info = exp.A[i]        
        samps.append('A' + str(i+1))
        conds.append('0')
        cmds.append(output + '/A' + str(i+1) + ' ' + info['First'] + ' ' + info['Second'])

    # For each replicate for mixture B...
    for i in range(len(exp.B)):
        samps.append('B' + str(i+1))
        conds.append('1')
        cmds.append(output + '/B' + str(i+1) + ' ' + info['First'] + ' ' + info['Second'])

    p = ThreadPool(6) 
    p.map(runKallisto, cmds)
    
    print 'Kallisto quantification completed'
    
    #
    # Reference: https://www.biostars.org/p/143458/#148465
    #            https://www.biostars.org/p/157240/#157242
    #
    
    # Do we have at least two replicates for each mixture?
    if len(exp.A) > 1 and len(exp.B) > 1:
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

    print 'RnaQuin quantification completed'

#
# VarQuin currently supports Kallisto and Salmon.    
#
def VarKReport(file):
    pass

class Experiment:
    def __init__(self, file):
        data = pandas.read_csv(file, sep='\t')        
        data['Mixture'] = data['Mixture'].astype('category')
        
        # Eg: 'A','A','B','B'
        mix = data['Mixture']
        
        if (len(mix.unique()) != 1 and len(mix.unique()) != 2):
            raise Exception('Invalid mixture column: ' + mix)        
        elif (not mix.isin(['A', 'B']).all()):
            raise Exception('Invalid mixture column. The values can only be A or B.')

        self.A = []
        self.B = []

        for index, row in data.iterrows():
            if row[0] == 'A':
                self.A.append({'First':row[1], 'Second':row[2]})
            else:
                self.B.append({'First':row[1], 'Second':row[2]})

#
# This file is part of Anaquin and not designed for external usage. Error checking is kept to minmial.
#
#   Eg: python ~/Sources/QA/scripts/kReport.py RnaQuin ARN024_v001.index /tmp/kallisto experiment.txt
#   Eg: python kReport.py RnaReport RnaKReport_report.pdf /tmp/kallisto
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
        RnaQuin(file, output, Experiment(sys.argv[4]))

    elif __mode__ == 'VarQuin':
        VarQuin(file, output, pair1, pair2)

    elif __mode__ == 'RnaReport':        
        createRNReport(file, output)

    elif __mode__ == 'VarReport':
        pass