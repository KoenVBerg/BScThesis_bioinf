#! /usr/bin/env python


"""
Author: Koen van den Berg
Function: This script performs multipe allignment using clustalW and calculates
the distance matrix using distmat. The Jukes-Cantor method is used for
correction as this is the most general method for allignment.  

Input:
(1) A file containing fasta sequences that are ready to be alligned.

Output:
(1) An .aln file, which contains the allingment file from clustalw
(2) An .distmat file, which contains the distance matrix of the allignment 
"""

# Import statements:
import os.path
import subprocess
from sys import argv
import sys
import argparse

# Functions:
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="This script performs a \
    multiple allignment and outputs a distance matrix.")
    parser.add_argument("-f", "--infiles", nargs='+',help="Name the input \
    files here", required = True)
    return(parser.parse_args())

def run_allignment(Infile):
    """
    Runs the allignment using clustalW
    Infile -- string, name of the user input file, fasta format
    """
    cmd_allignment = "clustalw -infile={}".format(Infile)
    cmd_rmtmp = "rm {}.dnd".format(Infile[:-6])
    res_aln = subprocess.check_output(cmd_allignment, shell=True)
    res_rm = subprocess.check_output(cmd_rmtmp, shell=True)

    return()    


def distance_matrix(alnfile):
    """
    Makes the distance matrix from the allignment file using distmat
    alnfile -- string, name of the allignment file, .aln
    """
    cmd_distmat = "distmat -sequence {} -nucmethod 1 -outfile \
    {}.distmat".format(alnfile, alnfile[:-4])
    res_distmat = subprocess.check_output(cmd_distmat, shell=True)
    return()


if __name__ == "__main__":

    # Obtaining the user arguments
    User_input = get_arguments()
    
    # Run the allignment:
    for file in User_input.infiles:
        print("Runs the allignment for {}".format(file))
        run_allignment(file) 
        Allignment_file = "{}.aln".format(file[:-6])
        distance_matrix(Allignment_file)



