#
# This file implements the K-mer expression pipeline for Anaquin. Users are expected to install the dependencies and thus not
# supported by Anaquin.
#
#   - Check installation
#   - Execute the pipeline
#

import os
import sys

def run(cmd):
    print(cmd)
    os.system(cmd)

def which(program):
    import os
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

def RNAKExpress(index, output, pair1, pair2):
    if which('./kallisto') is None:
        raise Exception('Kallsito not installed. Please check and try again.')

    quant = './kallisto quant -i ' + index + ' -o ' + output + ' ' + pair1 + ' ' + pair2

    # Quantify the input sequence files
    run(quant)
    
    print 'Kallisto completed'
    
def VarKExpress(file):
    pass

#
# Eg: python kexpress.py RnaQuin ARN004.v032.index /tmp/kallisto LRN087.1_val_1.fq LRN087.2_val_2.fq
#
if __name__ == '__main__':

    # Which sequin?
    mode = sys.argv[1]

    # Eg: A.R.3.index
    index = sys.argv[2]

    # Where the outputs should be written to
    output = sys.argv[3]

    # Sequence file for the first pair
    pair1 = sys.argv[4]
    
    # Sequence file for the second pair
    pair2 = sys.argv[5]
    
    print (mode)
    print (index)
    print (output)
    print (pair1)
    print (pair2)

    if not os.path.exists(output):
        os.makedirs(output)    

    if mode == 'RnaQuin':
        RNAKExpress(index, output, pair1, pair2)
    elif mode == 'VarQuin':
        VarKExpress(index, output, pair1, pair2)
