#
# This file implements the k-mer quantification pipeline for Anaquin.
#

import io
import os
import sys
import tempfile
import subprocess

def run(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        return e.output

def createText(file, txt):
    with io.FileIO(file, "w") as file:
        file.write(txt)
        
def runRScript(file):
    return run('Rscript ' + file)

def checkExec(program):
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
    
def RPackInstalled(pack):
    tmp = tempfile.gettempdir() + '/RPackInstalled.R'    
    cmd = pack + ' %in% rownames(installed.packages())'

    # Create an R-script for checking the existance of the package
    run('echo "require(' + pack + ')" > ' + tmp)

    # Eg: object 'ggplot2s' not found
    out = run('Rscript ' + tmp)
    
    return not ('no package' in out)

def RNAKReport(index, output, A1, A2, B1, B2):
    if checkExec('kallisto') is None:
        raise Exception('Kallisto is not installed. Please consult the user guide on www.sequin.xyz and try again.')

    if checkExec('R') is None:
        raise Exception('R is not installed. Please consult the user guide on www.sequin.xyz and try again.')

    if not RPackInstalled('sleuth'):
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

    # Do we have at least two replicates in each sample?
    #if len(A1) > 0 and len(B1) > 0:
    #    sleuth = """
    #             library(sleuth)
    #
                 # Where the Kallisto files are
    #             path <- '%s'
    #
                 # List of sample IDs
    #             sample_id <- dir(file.path(path))
    #             
                 # List of paths to the Kallisto results indexed by the sample IDs
    #             kal_dirs <- sapply(sample_id, function(id) file.path('/tmp/kallisto', "results", id, "kallisto"))
    #             
    #             samps <- c(%s)
    #             conds <- c(%s)
    #
                 # Construct full path for the samples
    #             path <- paste(path, samps, sep='/')
    #
    #             s2c <- data.frame(sample=samps, condition=conds)
    #             s2c <- dplyr::mutate(s2c, path=path)
    #
    #             so <- sleuth_prep(s2c, ~condition)
    #             so <- sleuth_fit(so)
    #             so <- sleuth_wt(so, 'conditionB')
    #
    #             results <- sleuth_results(so, 'conditionB')
    #
    #             write.csv(results, file='%s', row.names=FALSE, quote=FALSE)"""
    #             
    #    samps = ["'" + s + "'" for s in samps]
    #    conds = ["'" + s + "'" for s in conds]        
    #
    #    sleuth = sleuth % (output, ','.join(samps), ','.join(conds), "1")
    #
        # Create a R script for Sleuth
    #    createText(output + '/sleuth.R', sleuth)
    #
        # Run Sleuth
    #   runRScript(output + '/sleuth.R')

def VarKReport(file):
    pass

#
# Eg: python kReport.py RnaQuin ARN024.v032.index /tmp/kallisto 1,0 LRN087.1_val_1.fq LRN087.2_val_2.fq
#
if __name__ == '__main__':

    # Which sequin?
    mode = sys.argv[1]

    # Eg: A.R.3.index
    index = sys.argv[2]

    # Where the outputs should be written to
    output = sys.argv[3]
    
    if not os.path.exists(output):
        os.makedirs(output)    

    if mode == 'RnaQuin':
        
        # What's the design for A? Eg: 3,3
        A = 2*int(sys.argv[4].split(',')[0]) # Paired-end
    
        # What's the design for B? Eg: 3,3
        B = 2*int(sys.argv[4].split(',')[1]) # Paired-end

        A1 = []
        A2 = []
        B1 = []
        B2 = []
    
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

        RNAKReport(index, output, A1, A2, B1, B2)
    elif mode == 'VarQuin':
        VarKReport(index, output, pair1, pair2)
