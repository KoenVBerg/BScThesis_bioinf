#! usr/bin/env python3

"""
Author: Koen van den Berg
Function: This script takes the genbank files as input and makes a 
fasta file out of it. 

Explanation:
This script is designed to parse the genbank files of the pathways and to
produce fasta files of them. 

Input: 
(1) The genbank files, can be a space seperated list.
(2) Title, the name of the output file

Output:
(1) Output file which contains the fasta sequences
"""

# Import statements:
import os.path
import subprocess
from sys import argv
import argparse

# Functions:
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="This script takes the \
    genbank files as input and makes a fasta file out of it.")
    parser.add_argument("-f", "--infiles", help="Specify the names of the \
    input files here. This can be a space seperated list.", required = True, nargs = "+")
    parser.add_argument("-t", "--title", help="The name of the output file. \
    default = Pathways.fasta", required = True) 
    return(parser.parse_args())

def parse_genbank(Input_file):
    """
    Parses the genbank files specified by the user
    Input_file -- string, name of the input gb file
    """    
    Locus, Organism, DNA = "", "", []
    Flag = False
    with open(Input_file, "r") as f:
        for line in f:
            line = line.strip()
            if "LOCUS" in line:
                Locus = line.split()[1]
            if "SOURCE" in line:
                Organism = " ".join(line.split()[1:])
            if "ORIGIN" in line:
                Flag = True
            if Flag:
                DNA.append("".join(line.split()[1:]))
            if "//" in line:
                Flag = False
        yield(Locus, Organism, "".join(DNA))

if __name__ == "__main__":
    
    # Define user input:
    User_input = get_arguments()
    
    # Parse the genbank file
    Locus_l, Organism_l, DNA_l = [], [], []
    for File in User_input.infiles: 
        for Locus, Organism, DNA in parse_genbank(File):
            Locus_l.append(Locus)
            Organism_l.append(Organism)
            DNA_l.append(DNA)
    
    # Writing the Fasta file:
    Outfile = "Pathways.fasta"
    if User_input.title:
        Outfile = User_input.title
    with open(Outfile, "w") as o:
        for i in range(len(Locus_l)):
            o.write(">{}|{}\n".format(Locus_l[i], Organism_l[i]))
            o.write("{}\n".format(DNA_l[i]))
    
