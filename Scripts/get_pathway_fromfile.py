#! usr/bin/env python3

"""
Author: Koen van den Berg
Function: This script takes the output of the multigeneblast results as 
a .txt file and downloads the fasta files of the needed pathway results 
(Homologs). 

Explanation:
This script is designed to obtain the pathways from the multigeneblast results
quickly. Looking through the .html page of the results it seems rather hard to
download the pathway sequence from the html page directly, but this script
manages to do so in an easy manner. 

In summary, this script uses the multigeneblast results output.txt to obtain the
coordinates of the pathway. These are then used to parse out the sequence from
the whole genome sequence and parsed to fasta format.

First the multigeneblast result file, which is in text format, gets parsed to
obtain the accesion, organism and gene-list. The gene list is in a arbitrary
order. This means that to obtain the proper start and ending coordinates the
genelist should be sorted. This is being done with the use of lambda, which
enables sorting by the last 4 digits of the gene name. This results in an
ordered gene list. 

Next, the start and end coordinates from the first and last gene in the list are
obtained. This is done by using the efetch.entrez module to obtain the handle
(which is like opening a file) which contains the genbank file of the gene. In
this genbank file the coordinates are found and are parsed out. The coordinates
for both the first and the last gene in the sorted gene list are used because it
is not clear yet if the pathway is on the forward or backward DNA strand. 

Then, using the coordinates the pathway sequence is parsed out from the whole
genome sequence of the accession. It is checked whether the pathway is present
on the forward or backward DNA strand by checking which coordinate comes first.
Then, just like the genbank file, the Entrez module is used to obtain the wgs
sequence of the accession. 

At last a fasta file is written from the pathway sequence. 

Input:
(1) Number from the multigeneblast output of which you want the fasta sequence.
(2) Title of the pathway, this will be used in the title of the fasta sequence.
(3) Output file of the multigeneblast results as input for this script.

Output:
(1) Verbose report during the running of the script
(2) The fasta sequence of the pathway selected. 
"""

# Import statements:
import argparse
from Bio import Entrez, Seq, SeqIO

# Functions:
def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="Downloads the \
    pahtway files of the accession specified using the .txt output \
    file of the multigeneblast results.")
    parser.add_argument("-n", "--number", help="integer specifying the\
    number of the desired output on the results page.", required = True)
    parser.add_argument("-t", "--title", help="Put the title of the \
    pathway here", required = True)
    parser.add_argument("-i", "--input_file", help="The name of the\
    output file of the multigeneblast results that is used here as\
    input", required = True)
    return(parser.parse_args())

def parse_gene(File, number):
    """
    Parses multigeneblast output file and returns accession and 
    genelist
    File -- string, name of the multigeneblast output file
    number -- int, output number
    """
    Flag = False
    Hit = False
    Table = False
    List = []
    with open(File, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Details:"):
                Flag = True
            if Flag == True and line.startswith("{}.".format(number)):
                Hit = True
                accession = line.split()[1][:-2]
            if Hit == True and line.startswith("Source:"):
                s,x,contents = line.partition(":")
                name = contents.split(",")[0]
                organism = name.replace(" ", "_")
            if Hit == True and line.startswith("Table of Blast hits"):
                Table = True
            if Table == True and "\t" in line:
                Gene = line.split()[1]
                List.append(Gene)
            if Table == True and ">>" in line:
                Table = False
                Hit = False
    return(accession, organism, List)

def get_genbank(acc):
    """
    Obtains the genbank file for the accession and writes to tmp file
    acc -- string, the name of the accession
    """
    Entrez.email = "koen.vandenberg@wur.nl"
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="txt")
    tmp_file = "tmp.txt"
    with open(tmp_file, "w") as tmp:
        tmp.write(handle.read())
    handle.close()
    return(tmp_file)

def parse_coords(tmp_file):
    """
    parses the tmp file and obtains coords for the gene
    tmp_file -- string, the name of the tmp file
    """
    with open(tmp_file, "r") as tmp:
        for line in tmp:
            line = line.strip()
            if "/coded_by" in line:
                if "complement" in line:
                    junk,x,coords = line.partition(":")
                    start_coord = coords.split("..")[0]
                    end_coord = coords.split("..")[1][:-2]
                    return(int(start_coord), int(end_coord))
                else:
                    junk,x,coords = line.partition(":")
                    start_coord = coords.split("..")[0]
                    end_coord = coords.split("..")[1][:-1]
                    return(int(start_coord), int(end_coord))

def get_WGS(acc_code):
    """
    Downloads the whole genome sequence of the specified accession
    acc_code -- string, the name of the accession
    """
    Entrez.email = "koen.vandenberg@wur.nl"
    handle = Entrez.efetch(db="nucleotide", id=acc_code, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    #print(handle.readline().strip())
    handle.close()
    print(record)
    return(record.seq)


if __name__ == "__main__":
    
    # User input:
    User_input = get_arguments()
    
    # Obtaining the accession, organism and gene names:
    Accession, Organism, Gene_l = parse_gene(User_input.input_file, User_input.number)

    # Sorting the gene list 
    Gene_l_sorted = sorted(Gene_l, key = lambda x: x[-4:])
    
    # Acquiring the coordinates of the start and end gene:
    tmp_file = get_genbank(Gene_l_sorted[0])
    s_co_1, e_co_1 = parse_coords(tmp_file)
    
    tmp_file = get_genbank(Gene_l_sorted[-1])
    s_co_2, e_co_2 = parse_coords(tmp_file)
    print(s_co_1, e_co_1, s_co_2, e_co_2)
    
    # Get WGS and parse using coords
    WGS = get_WGS(Accession)
    
    if s_co_1 < e_co_2:
        Pathway_sequence = WGS[s_co_1:e_co_2]
    else:
        Pathway_sequence = WGS[s_co_2: e_co_1]

    # Writing the Fasta file:
    Title = User_input.title
    Outfile = "HOMOLOG_{}_{}_{}.fasta".format(Title, User_input.number, 
    Organism)
    with open(Outfile, "w") as o:
        o.write(">{}\n".format(Outfile[:-6]))
        o.write("{}\n".format(Pathway_sequence))

