#!/usr/bin/python

#
# This script generates a IGV session for sequin analysis
#

import os
import sys
import urllib

# URL for the in-silico chromosome
silicoFA = 'http://www.anaquin.org/downloads/chromo/1.0.0/chrT-2.1.0.tar.gz'

# Returns the URL for the transcriptome annotation
transGTF = 'http://www.anaquin.org/downloads/trans/TransStandard_1.0.gtf'

#
# Default template for generating a IGV session. We always show a in-silico chromosome, and it's assumed have a file name of chrT.fa.
#

sessionT = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{0}" hasGeneTrack="false" hasSequenceTrack="true" locus="chrT:1-44566700" path="session.xml" version="8">
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
def index(path, align):    
    print('Generating index for ' + align + ' to ' + path)
    return

    tmp = 'TEMP'

    # Silently create a temporary index directory
    os.system('mkdir -p ' + tmp)

    # Eg: /home/tedwong/ABCD.sam to ABCD.sam
    file = os.path.basename(align)

    # Eg: ABCD.sam to ABCD
    base = os.path.splitext(os.path.split(file)[-1])[0]

    # Create a sorted alignment, this is always needed for generating an index
    run('samtools sort ' + align + ' ' + tmp + '/' + base)

    # Generate the index
    run('samtools index ' + tmp + '/' + file)

    # Copy the alignment file
    run('cp ' + tmp + '/' + file + ' ' + path)

    # Move the index
    run('mv ' + tmp + '/*.bai ' + path)
    
    # Cleanup the index directory
    run('rm -rf ' + tmp)

def retrieve(url, path, file):
    return
    try:
        file = path + '/' + file
        
        #print('Downloading ' + url + ' to ' + file)
        urllib.urlretrieve(url, file)
        
        if file.endswith('tar'):
            run('tar -xvf ' + file + ' -C ' + path)
        elif file.endswith('tar.gz'):
            run('tar -zxvf ' + file + ' -C ' + path)
            run('rm -rf ' + file)

    except e:
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

    # Update the primary input (usually an alignment file)
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

def generateFusion(path, align):
    files = [ silicoGTF,
              'http://www.anaquin.org/downloads/fusion/normal_genes.gtf',
              'http://www.anaquin.org/downloads/fusion/fusion_genes.gtf']

    index(path, align)
    session(path, align, download(path, files))

def generateLadder(path, align):
    index(path, align)
    session(path, align, download(path, []))

def generateVar(path, align):
    files = [ 'http://www.anaquin.org/downloads/variant/DNA.variant.bed' ]

    index(path, align)
    session(path, align, download(path, files))

def generateTrans(path, align):
    index(path, align)
    session(path, align, download(path, [ silicoFA, transGTF ]))

if __name__ == '__main__':

    # Where to generate    
    mode = sys.argv[1]

    path = sys.argv[2]    
    path = os.path.dirname(os.path.abspath(path)) + '/' + path

    # Silently create the directory
    os.system('mkdir -p ' + path)
        
    if (os.path.isdir(path) == False):
        raise Exception('Invalid directory: ' + path)
    
    if (mode == 'Fusion'):
        generateFusion(path, sys.argv[3])
    elif (mode == 'Trans'):
        generateTrans(path, sys.argv[3])
    elif (mode == 'Variant'):
        generateVar(path, sys.argv[3])
    elif (mode == 'Ladder'):
        generateVar(path, sys.argv[3])

