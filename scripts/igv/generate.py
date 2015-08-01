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

import os
import sys
import urllib

def generate(path):

    # Create a IGV folder
    os.system('mkdir -p ' + path)

    #
    # Download the required files
    #

    urllib.urlretrieve("http://www.anaquin.org/downloads/chromo/chrT.fa", path + "chrT.fa")
    urllib.urlretrieve("http://www.anaquin.org/downloads/fusion/fusion_genes.gtf", path + "fusion_genes.gtf")
    urllib.urlretrieve("http://www.anaquin.org/downloads/fusion/normal_genes.gtf", path + "normal_genes.gtf")
    urllib.urlretrieve("http://www.anaquin.org/downloads/transcriptome/chrT_rna.gtf", path + "chrT_rna.gtf")

    #
    # Generate a IGV session file for the files
    #

    with open ("scripts/igv/session.template", "r") as f:
        session = f.read()

    with open(path + "/igv_session.xml", "w") as f:
        f.write(session)

if __name__ == '__main__':
    generate('/Users/tedwong/Sources/QA/A/')