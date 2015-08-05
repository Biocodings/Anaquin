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
silicoFA = 'www.anaquin.org/downloads/chromo/chrT.fa'

# Returns the URL for the in-silico annotation
silicoGTF = 'www.anaquin.org/downloads/transcriptome/chrT_rna.gtf'

#
# Default template for generating a IGV session. We always show a in-silico chromosome, and it's assumed have a file name of chrT.fa.
#

template = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="{0}/chrT.fa" hasGeneTrack="false" hasSequenceTrack="true" locus="chrT:1-44566700" path="/Users/tedwong/Sources/QA/A/igv_session.xml" version="8">
    <Resources>
        {1}
    </Resources>
    <Panel height="574" name="FeaturePanel" width="1423">
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/>
        <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="/Users/tedwong/Sources/QA/A/chrT_rna.gtf" name="chrT_rna.gtf" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>
    </Panel>
    <PanelLayout dividerFractions="0.010256410256410256"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>"""

# Generate an index for a SAM/BAM file
def index(path, align):

    # Eg: accepted_hits
    index = os.path.splitext(os.path.split(align)[-1])[0]
    
    # Generate the index
    os.system('samtools index ' + align + ' ' + path + "/" + index)

# Download the required files and generate a IGV session 
def download(path, files):

    # Create a session folder
    os.system('mkdir -p ' + path)

    #
    # Download the required files
    #

    #urllib.urlretrieve("http://www.anaquin.org/downloads/chromo/chrT.fa", path + "/chrT.fa")
    urllib.urlretrieve("http://www.anaquin.org/downloads/fusion/fusion_genes.gtf", path + "fusion_genes.gtf")
    urllib.urlretrieve("http://www.anaquin.org/downloads/fusion/normal_genes.gtf", path + "normal_genes.gtf")
    urllib.urlretrieve("http://www.anaquin.org/downloads/transcriptome/chrT_rna.gtf", path + "chrT_rna.gtf")

    return ["chrT_rna.gtf"]

def session(path, files):
    global template
    
    # IGV assumes a full path
    path = os.path.abspath(path) + '/'

    # Eg: <Resource path="{0}/chrT_rna.gtf"/>
    res = ''

    for file in files:        
        res = res + '<Resource path="' + path + file + '"/>\n'        

        #
        # Eg: <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="/Users/tedwong/Sources/QA/A/chrT_rna.gtf" name="chrT_rna.gtf" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>
        #

    # Update the specifed directory
    template = template.replace('{0}', path)

    # Update the specifed directory
    template = template.replace('{1}', res[:-1])

    with open(path + "/igv_session.xml", "w") as f:
        f.write(template)

def generateFusion(path, align):
    files = [ silicoFA,
              silicoGTF,
              'http://www.anaquin.org/downloads/fusion/normal_genes.gtf',
              'http://www.anaquin.org/downloads/fusion/fusion_genes.gtf']

    #index(path, align)
    session(path, download(path, files))

if __name__ == '__main__':

    # Where to generate    
    mode = sys.argv[1]

    # Where the files should be generated    
    path = sys.argv[2]
    
    if (mode == 'Fusion'):
        generateFusion(path, sys.argv[3])

