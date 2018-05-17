 #! usr/bin/env python3

"""
Author: Koen van den Berg
Function: This script checks all the sequences using fastqc. It uses 
the "fastqc_data.txt" file to check if the sequences in the sample are 
good enough for further usage in the pipeline. If this is not the case 
the file will be trimmed using trimmomatic. 

Explanation:
This script is designed to check the fastq files for poor quality sequences.
This is quite important, because in the end result those poor sequences should
not be added to the results. Thus, this script uses the command line program
fastqc to asses whether the fastq files are of appropriate quality. It does so
by checking all the inputted fastq-files and producing one big summary table of
the results. In this summary table the number of poorly flagged sequences is
mentioned. 

Input:
(1) fastq files, can be space seperated list 

Output:
(1) a summary table that contains the Filename, Nseq, Npoor, Seqlen and GC_content.

"""

# Import statements:
import argparse
import os.path
import subprocess
from sys import argv

# Functions:
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="This script checks all \
     the sequences using fastqc. It uses the fastqc_data.txt file to check \
     if the sequences in the sample are good enough for further usage in the \
     pipeline. If this is not the case the file will be trimmed using \
    trimmomatic. the output file of this script is Summary_fastqc.txt")
    parser.add_argument("-f", "--infiles", help="the input fastq-files that need \
    to be checked. This can be a space seperated list.", required = True, nargs= "+")
    return(parser.parse_args())
 

def fastqc_check(fastq_file):
    """
    runs the fastqc program on the accession and extracts the output 
    zip.
    fastq_file -- String, the name of the accesion sample, ".fastq"-file
    """
    Path = fastq_file[:-6] + "_fastqc/"
    if os.path.exists(Path):
        print("the fastqc quality check has already been performed! Proceeds \
to make summary table.")
        return()
    else:
        cmd_fastqc = "fastqc {} --extract".format(fastq_file)
        res_download = subprocess.check_output(cmd_fastqc, shell=True)
        return()

def parse_data(fastq_file):
    """
    Parses the fastqc output file 
    fastq_file -- String, the name of the accesion sample, ".fastq"-file
    """    
    File =  "{}_fastqc/fastqc_data.txt".format(fastq_file[:-6])
    filename, nseq, npoor, seqlen, GC = "", 0, 0, 0, 0
    with open(File, "r") as f:
        for line in f:
            line = line.strip()
            if "Filename" in line: 
                filename = line.split()[1]
            if "Total Sequences" in line: 
                nseq = line.split()[2]
            if "Sequences flagged as poor quality" in line:
                npoor = line.split()[5]
            if "Sequence length" in line:
                seqlen = line.split()[2]
            if "%GC" in line:
                GC = line.split()[1]
        yield(filename, nseq, npoor, seqlen, GC)

if __name__ == "__main__":
    
    # Define user input:
    User_input = get_arguments() 
    
    # Perform fastqc check
    for fqfile in User_input.infiles:
        fastqc_check(fqfile)
        
    # Parse summary files to list:
    finame_l, nseq_l, npoor_l, seqlen_l, GC_l = [], [], [], [], []
    for fqfile in User_input.infiles:
        for filename, nseq, npoor, seqlen, GC in parse_data(fqfile):
            finame_l.append(filename)
            nseq_l.append(nseq)
            npoor_l.append(npoor)
            seqlen_l.append(seqlen)
            GC_l.append(GC)
    
    # Write lists to summary table:
    Outfile = "Summary_fastqc.txt"
    with open(Outfile, "w") as f:
        f.write("Filename\tNseq\tNpoor\tSeqlen\tGC_content\n")
        for i in range(len(finame_l)):
            f.write("{}\t{}\t{}\t{}\t{}\n".format(finame_l[i], nseq_l[i], 
            npoor_l[i], seqlen_l[i], GC_l[i]))
        
