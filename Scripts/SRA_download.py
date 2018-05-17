#! usr/bin/env python3

"""
Author: Koen van den Berg
Function: This script downloads all the sra files from the NCBI online 
sample database. Next, if specified, the fastq files can be extracted 
from the obtained .sra files using fastqdump. 

Explanation script: 
This script is designed to download and acquire sra files from the ncbi 
website. On this site the data used by most research papers is available 
and can be acquired using very simple python syntax. In this script the 
sra files, which are the packaged file, are downloaded fromt this database. 
The link to the website on the ncbi where the files can be found is:
ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP 

Using this link the sra files can be downloaded with wget. Next, the fastq 
files are extracted from the sra files using the program fastq-dump. This 
program can be accessed from the command line (like usual) and is used in this
script to extract the fastq files from the sra files.

Input: 
(1) a .txt file which contains the run(s) found on the ncbi website. This .txt
file should contain 1 run per line.
(2) The name of the study as found on the ncbi website of that paper, for 
example: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA188481&go=go. 

Output:
(1) The downloaded SRA files and 
(2) if enabled the fastq files which are extracted from the .sra files.   
"""

# Import statements
import os.path
import subprocess
from sys import argv
import argparse

# Functions
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="This script \
    downloads all the sra files from the NCBI online sample database. \
    Next, if specified, the fastq files can be extracted from the \
    obtained .sra files using fastqdump.")
    parser.add_argument("-i", "--infile", help="File containing the sra \
    accession numbers that you want to download. Should contain one \
    sra code per line", required = True)
    parser.add_argument("-s", "--study", help="The study code from \
    which you want to download the sra files. This code can be found at \
    ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP \
    ", required = True)
    parser.add_argument("-f", "--fastq", help="Boolean variable \
    specifying if the fastq files should be extracted from the \
    downloaded .sra files. Options: True/False", default = False)
    return(parser.parse_args())
    
    
def get_acc(File):
    """
    Yields the Accesion numbers from the user specified file
    file -- string, name of the input file
    """
    with open(File, "r") as f:
        for line in f:
            line = line.strip()
            yield(line)
        
def download_sra(Accession, study_code):
    """
    Downloads the RSA file of the specified the Accesion number
    Accession -- String, the name of the accesion sample
    study_code -- string, the codename of the study
    """
    SRA_file = str(Accession + ".sra")
    if os.path.exists(SRA_file):
        print("The SRA file already exists")
        return()
    else:
        Path = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/" + study_code[0:6] + "/" + study_code + "/" + Accession + "/" + SRA_file
        cmd_download = "wget" + " " + Path
        res_download = subprocess.check_output(cmd_download, shell=True)
        return()

def make_fastq(Accession):
    """
    Makes a fastq file of the rsa file
    Accession -- String, the name of the accesion sample
    """
    SRA_file = str(Accession + ".sra")
    Fastq_file = str(Accession + ".fastq")
    if os.path.exists(Fastq_file):
        print("The fastq file already exists")
        return()
    else:
        cmd_fastq = "fastq-dump" + " " + SRA_file
        res_fastq = subprocess.check_output(cmd_fastq, shell=True)
        return()

if __name__ == "__main__":
    
    # Obtaining the user input arguments
    User_input = get_arguments()
   
    # Download the files and extract fastq-files  
    for Acc in get_acc(User_input.infile):
        download_sra(Acc, User_input.study)
        if User_input.fastq == True:
            make_fastq(Acc)
    
