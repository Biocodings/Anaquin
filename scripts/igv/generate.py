#!/usr/bin/python
#
# This script generates a IGV session for fusion analysis. Rather than keeping the required files in the distribution,
# the script will download the files online. This saves the bunden of distributing the in-silico chromosome.
#
#   - In silico chromosome
#   - Reference GTF for fusion genes
#   - Reference GTF for normal genes
#   - Reference GTF for the RNA standard
#
#   - Generated alignment file
#

import os
import sys
import urllib

# URL for the in-silico chromosome
silicoFA = 'www.anaquin.org/downloads/chromo/chrT.fa'

# Returns the URL for the in-silico annotation
silicoGTF = 'www.anaquin.org/downloads/transcriptome/chrT_rna.gtf'

# Generate an index for a SAM/BAM file
def generateIndex(path, align):
    # Generate the index
    os.system('samtools index ' + align + ' ' + os.path.splitext(align)[0])

    # Move the index files to the specified path (samtools creates index to the directory for the alignment file)
    os.system('mv index.pbi ' + path)

# Download the required files and generate a IGV session 
def generate(path, files):

    # Create a IGV folder
    os.system('mkdir -p ' + path)

    #
    # Download the required files
    #

    urllib.urlretrieve("http://www.anaquin.org/downloads/chromo/chrT.fa", path + "/chrT.fa")
    urllib.urlretrieve("http://www.anaquin.org/downloads/fusion/fusion_genes.gtf", path + "/fusion_genes.gtf")
    urllib.urlretrieve("http://www.anaquin.org/downloads/fusion/normal_genes.gtf", path + "/normal_genes.gtf")
    urllib.urlretrieve("http://www.anaquin.org/downloads/transcriptome/chrT_rna.gtf", path + "/chrT_rna.gtf")

    #
    # Generate a IGV session file for the files
    #

    with open ("session.template", "r") as f:
        session = f.read()

    with open(path + "/igv_session.xml", "w") as f:
        f.write(session)

def generateFusion(path, align):
    files = [ silicoFA,
              silicoGTF,
              'http://www.anaquin.org/downloads/fusion/normal_genes.gtf',
              'http://www.anaquin.org/downloads/fusion/fusion_genes.gtf']

    print('Download Files')
    #generate(path, files)
    
    print('Generate Index')
    generateIndex(path, align)

if __name__ == '__main__':

    # What to generate    
    mode = sys.argv[1]

    # Where the files should be generated    
    path = sys.argv[2]

    if (mode == 'Fusion'):
        generateFusion(path, sys.argv[3])



