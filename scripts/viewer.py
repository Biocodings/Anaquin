#!/usr/bin/python

#
# This script generates an IGV session for sequin analysis. It is a part of the Anaquin distribution, where it's being bundled with
# the C++ code. The script can also run by itself.
#

import os
import sys
import urllib

#
# TransQuin
#

T_CHR_T = 'https://s3.amazonaws.com/anaquin/chromosomes/CTR001.v021.fa'
T_STAND = 'https://s3.amazonaws.com/anaquin/annotations/ATR001.v032.gtf'

#
# VarQuin
#

# URL for the in-silico chromosome (variant)
V_CHR_T  = 'https://s3.amazonaws.com/anaquin/chromosomes/CVA002.v015.fa'

# URL of the variant standard
V_STAND = 'https://s3.amazonaws.com/anaquin/annotations/AVA017.v032.bed'

# URL of the known variants
V_VAR = 'https://s3.amazonaws.com/anaquin/annotations/AVA009.v032.vcf'

# URL of the normal standard
fusNStand = 'https://s3.amazonaws.com/anaquin/annotations/????'

# URL of the fusion standard
fusFStand = 'https://s3.amazonaws.com/anaquin/annotations/????'

# URL of the synthetic metagenomic community
commAmazon = 'https://s3.amazonaws.com/anaquin/chromosomes/CME003.v013.fa'

# URL of the metagenomic standard
mStandAmazon = 'https://s3.amazonaws.com/anaquin/annotations/AME015.v032.bed'

#
# Default template for generating a IGV session. We always show a in-silico chromosome, and it's assumed have a file name of chrT.fa.
#

sessionT = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{0}" hasGeneTrack="false" hasSequenceTrack="true" path="session.xml" version="8">
    <Resources>
        {1}
        {2}
    </Resources>
    <Panel height="239" name="Panel1438759846196" width="1423">
        {ALIGN}
    </Panel>
    <Panel height="574" name="FeaturePanel" width="1423">
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/>
        {3}
    </Panel>
    <PanelLayout dividerFractions="0.010256410256410256,0.42735042735042733"/>    
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>"""

alignTrackT = """
        <Track altColor="0,0,178" autoScale="true" color="175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{PATH}/{FILE}_coverage" name="{FILE} Coverage" showReference="false" snpThreshold="0.2" sortable="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="60.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="{FILE}" name="{FILE}" showSpliceJunctions="false" sortable="true" visible="true">
            <RenderOptions colorByTag="" colorOption="UNEXPECTED_PAIR" flagUnmappedPairs="false" groupByTag="" maxInsertSize="1000" minInsertSize="50" shadeBasesOption="QUALITY" shadeCenters="true" showAllBases="false" sortByTag=""/>
        </Track>
"""

# Template for a custom track
trackT = """<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{0}{1}" name="{1}" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/> """

# Template for a resource such as GTF
resourceT = """        <Resource path="{1}"/> """

def run(cmd):
    os.system(cmd)

# Generate an index for a SAM/BAM file
def index(path, file):
    print('Generating index for ' + file + ' to ' + path)

    tmp = '/tmp'

    # Silently create a temporary directory for indexing
    os.system('mkdir -p ' + tmp)

    # Eg: /home/tedwong/ABCD.bam to ABCD.bam
    file = os.path.basename(file)

    # Eg: ABCD.sam to ABCD
    base = os.path.splitext(os.path.split(file)[-1])[0]

    if (file.endswith('.sam')):
        bam = file.replace('.sam', '.bam')
        
        run('samtools view -Sb ' + file + ' | samtools sort - ' + tmp + os.sep + base)
        run('samtools index ' + tmp + os.sep + bam)
        
        # From now on, we'll deal only with BAM format
        file = bam
    else:
        run('samtools sort ' + file + ' ' + tmp + '/' + base)
        run('samtools index ' + tmp + '/' + file)

    # Copy the alignment file
    run('cp ' + tmp + '/' + file + ' ' + path)

    # Move the index
    run('mv ' + tmp + '/*.bai ' + path)
    
    return path + '/' + file

def retrieve(url, path, file):
    try:
        file = path + '/' + file
        
        print('Downloading ' + url + ' to ' + file)
        urllib.urlretrieve(url, file)
        
        if file.endswith('tar'):
            run('tar -xvf ' + file + ' -C ' + path)
        elif file.endswith('tar.gz'):
            print('Extracting the file')
            run('tar -zxvf ' + file + ' -C ' + path)
            run('rm -rf ' + file)
    except:
        raise Exception('Failed to retrieve: ' + path)

# Download the required files and generate a IGV session 
def download(path, files):
    # Create a session folder
    os.system('mkdir -p ' + path)

    #
    # Download the required files
    #
    
    downloads = []

    for file in files:
        if (file.startswith('http')):
            retrieve(file, path, file.split('/')[-1])
        downloads.append(file.split('/')[-1])

    return downloads

def session(path, align, files):
    global resourceT, sessionT, trackT, alignTrackT

    # IGV assumes a full path
    path = os.path.abspath(path) + '/'

    res = ''    
    tra = ''

    for i in range(0, len(files)):
        if i>0:        
            tra  = tra + trackT.replace('{0}', path).replace('{1}', files[i])
            res  = res + resourceT.replace('{0}', path).replace('{1}', files[i]) + '\n'       

    #tra  = tra + trackT.replace('{0}', path).replace('{1}', align)
    res  = res + resourceT.replace('{0}', path).replace('{1}', align) + '\n'

    # We can assume the first file is always the reference genome and it's released in the GZ format
    chrT = path + os.path.basename(files[0]).replace('.tar.gz', '.fa')

    # Update the reference genome
    sessionT = sessionT.replace('{0}', chrT)

    # Update the alignment file)
    sessionT = sessionT.replace('{1}', '')
    
    #
    # Create a custom track with coverage for the alignment
    #

    alignTrackT = alignTrackT.replace('{FILE}', align)
    alignTrackT = alignTrackT.replace('{PATH}', path)
    sessionT    = sessionT.replace('{ALIGN}'  , alignTrackT)

    #
    # Create resoures for the custom tracks
    #

    sessionT = sessionT.replace('{2}', res[:-1])
    sessionT = sessionT.replace('{3}', tra[:-1])

    with open(path + "/session.xml", "w") as f:
        f.write(sessionT)

#
# Generate an IGV session for transcriptome analysis
#

def generateTrans(path, align, chrT=T_CHR_T, stand=T_STAND):
    
    align = index(path, align)
    session(path, align, download(path, [ chrT, stand ]))

#
# Generate an IGV session for variant analysis
#

def generateVar(path, align, chrT=V_CHR_T, stand=V_STAND, variant=V_VAR):
    align = index(path, align)
    session(path, align, download(path, [ chrT, stand, variant ]))

#
# Generate an IGV session for fusion analysis
#

def generateFusion(path, align, chrT=None, stand=None):
    align = index(path, align)
    session(path, align, download(path, [ chrTFA, transStand, fusNStand, fusFStand ]))

#
# Generate an IGV session for metagenomic analysis
#

def generateMeta(path, align):
    align = index(path, align)
    session(path, align, download(path, [ commAmazon, mStandAmazon ]))

def printUsage():
        print '\nProgram: viewer.py (Tool for generating IGV session)'
        print 'Version: 1.0.0\n'
        print 'Usage: python viewer.py <TransQuin|FusQuin|VarQuin|MetaQuin|LadQuin> <Output> <Alignment File>'
        print ''

if __name__ == '__main__':

    if (len(sys.argv) != 4):
        printUsage()
        exit()

    mode = sys.argv[1]
    path = sys.argv[2]
    file = sys.argv[3]
    
    if (not os.path.isfile(file)):
        raise Exception('Invalid file: ' + file)
    
    if (not os.path.isabs(path)):
        path = os.path.dirname(os.path.abspath(path)) + '/' + path    
        
    # Silently create the directory
    os.system('mkdir -p ' + path)

    if (os.path.isdir(path) == False):
        raise Exception('Invalid directory: ' + path)
    
    # Eg: accepted_hits.bam
    align = sys.argv[3]

    if (mode == 'TransQuin'):
        generateTrans(path, file)
    elif (mode == 'FusQuin'):
        generateFusion(path, file)
    elif (mode == 'VarQuin'):
        generateVar(path, file)
    elif (mode == 'MetQuin'):
        generateMeta(path, file)
    else:
        printUsage()