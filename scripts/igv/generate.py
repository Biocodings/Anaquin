#!/usr/bin/python
#
# This script generates a IGV session for fusion analysis. Rather than keeping the required files in the distribution,
# the script will download the files online. This saves the bunden of distributing the in-silico chromosome.
#
# Requirments:
#
#   - SAMTools
#   - Active Internet connection
#

import os
import sys
import urllib

# URL for the in-silico chromosome
silicoFA = 'http://www.anaquin.org/downloads/chromo/chrT.fa'

# Returns the URL for the in-silico annotation
silicoGTF = 'http://www.anaquin.org/downloads/transcriptome/RNA.v1.gtf'

#
# Default template for generating a IGV session. We always show a in-silico chromosome, and it's assumed have a file name of chrT.fa.
#

sessionT = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{0}chrT.fa" hasGeneTrack="false" hasSequenceTrack="true" locus="chrT:1-44566700" path="/Users/tedwong/Sources/QA/A/igv_session.xml" version="8">
    <Resources>
        <Resource path="accepted_hits.bam"/>        
        {1}
    </Resources>
    <Panel height="239" name="Panel1438759846196" width="1423">
        <Track altColor="0,0,178" autoScale="true" color="175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="/Users/tedwong/Sources/QA/scripts/Temp/accepted_hits.bam_coverage" name="accepted_hits.bam Coverage" showReference="false" snpThreshold="0.2" sortable="true" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="60.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="/Users/tedwong/Sources/QA/scripts/Temp/accepted_hits.bam" name="accepted_hits.bam" showSpliceJunctions="false" sortable="true" visible="true">
            <RenderOptions colorByTag="" colorOption="UNEXPECTED_PAIR" flagUnmappedPairs="false" groupByTag="" maxInsertSize="1000" minInsertSize="50" shadeBasesOption="QUALITY" shadeCenters="true" showAllBases="false" sortByTag=""/>
        </Track>
    </Panel>
    <Panel height="574" name="FeaturePanel" width="1423">
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/>
        {2}
    </Panel>
    <PanelLayout dividerFractions="0.010256410256410256,0.42735042735042733"/>    
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>"""

trackT = """<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="{0}{1}" name="{1}" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/> """

# Template for a resource such as GTF
resourceT = """        <Resource path="{1}"/> """

# Generate an index for a SAM/BAM file
def index(path, align):

    # Eg: accepted_hits
    index = os.path.splitext(os.path.split(align)[-1])[0]
    
    # Generate the index
    os.system('samtools index ' + align)

    # Copy the alignment file
    os.system('cp ' + align + ' ' + path)

    # Copy the alignment file
    os.system('cp ' + align + '.bai ' + path)

# Download the required files and generate a IGV session 
def download(path, files):

    # Create a session folder
    os.system('mkdir -p ' + path)

    #
    # Download the required files
    #
    
    downloads = []

    # We always need the chromosome
    urllib.urlretrieve(silicoFA, path + 'chrT.fa')    

    for file in files:
        urllib.urlretrieve(file, path + file.split('/')[-1])
        downloads.append(file.split('/')[-1])

    return downloads

def session(path, files):
    global resourceT, sessionT, trackT

    # IGV assumes a full path
    path = os.path.abspath(path) + '/'

    res = ''    
    tra = ''

    for file in files:
        
        #
        # {0}: specified directory
        # {1}: resource file
        #

        res  = res + resourceT.replace('{0}', path).replace('{1}', file) + '\n'       
        tra  = tra + trackT.replace('{0}', path).replace('{1}', file)        

    # Update the specifed directory
    sessionT = sessionT.replace('{0}', path)

    # Update the specifed directory
    sessionT = sessionT.replace('{1}', res[:-1])
    sessionT = sessionT.replace('{2}', tra[:-1])

    with open(path + "/igv_session.xml", "w") as f:
        f.write(sessionT)

def generateFusion(path, align):
    files = [ silicoGTF,
              'http://www.anaquin.org/downloads/fusion/normal_genes.gtf',
              'http://www.anaquin.org/downloads/fusion/fusion_genes.gtf']

    index(path, align)
    session(path, download(path, files))

if __name__ == '__main__':

    # Where to generate    
    mode = sys.argv[1]

    # Where the files should be generated    
    path = sys.argv[2]
    
    if (mode == 'Fusion'):
        generateFusion(path, sys.argv[3])

