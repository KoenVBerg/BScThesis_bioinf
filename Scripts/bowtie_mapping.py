#! usr/bin/env python3

"""
Author: Koen van den Berg
Function: This script maps the fastq files against the pathways (in 
fasta format). 
Usage: use -h flag

Explanation:
This script is designed to map the fastq files against a reference, in this case
the reference is the file containing all the pathways and bacterial reference
genes and the samples are the fastq-files. 

Overview steps:
(1) build reference index
(2) map query against reference index
(3) convert SAM --> BAM
(4) find pair and merge BAM if paired data
(5) sort BAM file
(6) index BAM file
(7) count BAM file
(8) parse count
(9) calculate TPM count
(10) output tabular file

First the bowtie2 index gets build (if not already build) using the
bowtie2-build command. The bowtie2 index gets stored in the same folded in which
the reference file is present. Hereafter the bowtie2 index will be used for the
mapping procedure. The mapping is done with local alignment enabled, as the
query sequences are each 200bp long. Eight processor units will be used during
the mapping step and a SAM file will be produced. This SAM file will then be
converted to a BAM file. 

If the query sequence is run on a paired lane in Illumia Hiseq, the -f flag
should be used in combination with the SRA runtable from the NCBI run selector
repository. This file will be used to find the pair of the fastq-file. Only if 
the BAM file is present the merging will be initiated. This will in the end
result in a merged BAM file which contains all the information for the paired
sample. Important to note here is that the current samples are only run on a 
seperate lane while not being actual paired end data. 

Next the BAM file will be indexed and counted using samtools index and samtools
idxstats. The resulting binary output string is stored into memory and
translated to unicode. Then the TPM counts will be calculated and outputted into
one big table. Take a look at the code below the if name == main statement.

Input: 
(1) reference FASTA file
(2) fastq query files 
(3) OPTIONAL: SraRuntable for paired data 
(4) name of the output file
(5) method of local alignment

Output: 
(1) SAM, BAM, sorted BAM and BAI files
(2) name_output_table.txt

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
    parser = argparse.ArgumentParser(description="This script maps the \
     fastq files to the reference pathways. After mapping the .sam file \
      gets converted to bam and sorted. Then the bam file gets indexed. \
      The samtools idxstats function is then able to count the reads \
      and outputs this in a binary string format. This string is then \
      parsed to counts which are normalized using the TPM formula. This \
      is then parsed to a neat count table that displays the samples and \
      pathways and corresponding counts.\
      Important!!!\
      Paired samples and unpaired samples cannot be run at the same \
      time")
    parser.add_argument("-r", "--reference", help="The fasta file that \
    includes the reference sequences. Should be in .fasta format", required = True)
    parser.add_argument("-i", "--fastq", nargs='+',help="This is the \
    input fastq sequences. This can be a space seperated list from \
    the command line", required = True)
    parser.add_argument("-f", "--paired", help="This should be a file \
    containing the metadata from the SRA site", default = False)
    parser.add_argument("-o", "--output", help="Specify the name of \
    the output file here", required = True)
    parser.add_argument("-m", "--method", help="very-fast-local; \
    fast-local; sensitive-local; very-sensitive-local", required = True)
    return(parser.parse_args())

def bowtie2_index(Ref):
    """
    Builds the index for the bowtie2 mapping
    Ref -- String, the name of the reference file, ".fasta"-file
    """
    Path = "{}.1.bt2".format(Ref[:-6])
    if os.path.exists(Path):
        print("Note: The Bowtie2 Index has already been created")
        return()
    elif not os.path.exists(Ref):
        print("Note: The Reference Fasta File does not exist")
        return()
    else:
        cmd_bowtie2_idx = "bowtie2-build {} {}".format(Ref, Ref[:-6])
        res_bwt2index = subprocess.check_output(cmd_bowtie2_idx, shell=True)
        return()

def bowtie2_map(Ref, Fq_file, Method):
    """
    Maps the .fq file to the reference (fasta)
    
    Ref -- String, the name of the reference file, ".fasta"-file
    Fq_file -- String, the name of the accesion sample, ".fastq"-file
    """
    Path = "{}.sam".format(Fq_file[:-6])
    if os.path.exists(Path):
        print("Note: The Bowtie2 mapping has already been performed \
        for this sample")
        return(Path)
    elif not os.path.exists(Fq_file):
        print("Note: The sample fastq-file does not exist")
        return()
    else:
        cmd_bowtie2_map = "bowtie2 -q -p 8 --{} -x {} -U {} -S \
        {}.sam".format(Method, Ref[:-6], Fq_file, Fq_file[:-6])
        res_map = subprocess.check_output(cmd_bowtie2_map, shell=True)
        return(Path)    
    
def make_bam(sam_file):
    """
    converts .sam to .bam using samtools view
    
    sam_file -- string, name of the outputted bowtie2 mapping
    """
    Path = "{}.bam".format(sam_file[:-4])
    if os.path.exists(Path):
        print("Note: The sam<->bam conversion has already been performed \
        for this sample")
        return(Path)
    elif not os.path.exists(sam_file):
        print("Note: The sample.sam file does not exist")
        return()
    else:
        cmd_bowtie2_bam = "samtools view -Sb {} > {}\
        ".format(sam_file, Path)
        res_bam = subprocess.check_output(cmd_bowtie2_bam, shell=True)
        return(Path)    

"""
def check_sample(Fq_file):
"""
    #Returns True if the fastq_file is paired
    
    #Fq_file -- String, the name of the accesion sample, ".fastq"-file
"""
    if int(Fq_file[-9]) < 5:
        return(True)
    if int(Fq_file[-9]) == 5:
        return(False)
"""    
    
    
def sort_bam(bam_file):
    """
    sorts the bam file
    
    bam_file -- String, the name of the accession bamfile, ".bam"-file
    """

    Path = "{}_sorted.bam".format(bam_file[:-4])
    if os.path.exists(Path):
        print("Note: The bam sorting has already been performed for this sample")
        return(Path)
    elif not os.path.exists(bam_file):
        print("Note: The sample.bam file does not exist")
        return()
    else:
        cmd_bam_sort = "samtools sort {} {}".format(bam_file, Path[:-4])
        res_sort = subprocess.check_output(cmd_bam_sort, shell=True)
        return(Path)  

        
def index_bam(sorted_bam_file):
    """
    Builds a bam index
    
    sorted_bam_file -- String, the name of the accesion sample, ".fastq"-file    
    """
    Path = "{}.bai".format(sorted_bam_file)
    if os.path.exists(Path):
        print("Note: The bam indexing has already been performed for this sample")
        return()
    else:
        cmd_bam_index = "samtools index {}".format(sorted_bam_file)
        res_index = subprocess.check_output(cmd_bam_index, shell=True)
        return()   
            
def count_bam(sorted_bam_file):
    """
    Counts the bam index file and returns binary string with counts and 
    pathways
    
    sorted_bam_file -- String, the name of sorted .bam file    
    """
    Path = sorted_bam_file
    cmd_bam_count = "samtools idxstats {}".format(Path)
    res_count = subprocess.check_output(cmd_bam_count, shell=True, 
    stderr = subprocess.STDOUT)
    return(res_count.decode("utf-8")) # Makes it a string i.s.o. binary
        
def parse_count(count_output):
    """
    Parses the count output string to a tabular format
    
    count_ouput -- String, the output of samtools    
    """
    contents = count_output.split("\n")
    for item in contents[:-2]:
        pathway = item.split("\t")[0]
        count = item.split("\t")[2]
        read_length = item.split("\t")[1]
        yield(read_length, pathway, count)
        
def calculate_rate_sum(lengths_l,counts_l):
    """
    Calculates the rate sum for a sample 
    
    lengths_l -- list, list of reference lengths
    counts_l -- list, list of counts to reference
    """    
    rate_sum = 0
    for i in range(len(lengths_l)): 
        rate = counts_l[i]/lengths_l[i]
        rate_sum += rate
    return(rate_sum)

def calculate_TPM(lengths_l,counts_l, rate_sum):
    """
    Calculates the rate sum for a sample 
    
    lengths_l -- list, list of reference lengths
    counts_l -- list, list of counts to reference
    """    
    TPM_counts = []
    for i in range(len(lengths_l)):
        rate = counts_l[i]/lengths_l[i]
        TPM = rate/rate_sum * 10**6
        TPM_counts.append(TPM)
    return(TPM_counts)

def find_pair(input_file, Fq_file):
    """
    Parses the SRA runtable and returns the sample pair
    
    Fq_file -- String, the name of the accesion sample, ".fastq"-file 
    input_file -- string, name of the SRA runtable file
    """
    # Finding the related experiment name
    Run, Experiment = "", ""
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            Run_1 = line.split("\t")[10]
            if Fq_file[:-6] == Run_1:
                Experiment_1 = line.split("\t")[2]
                break
    
    # Finding the paired sample:
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            Experiment_2 = line.split("\t")[2]
            if Experiment_1 == Experiment_2:
                Run_2 = line.split("\t")[10]
                if Run_1 != Run_2:
                    return(Run_2)
                

def merge_pair(fq_1, fq_2):
    """
    Merges the two paired bam files and returns merged bam file string
    
    fq_1 -- String, the name of the accesion sample, ".fastq"-file 
    fq_2 -- string, name of the accesion sample, does not have .fastq 
    extention!!!
    """
    
    Path_1 = "{}_{}.bam".format(fq_1[:-6], fq_2)
    Path_2 = "{}_{}.bam".format(fq_2, fq_1[:-6])
    bam_1 = "{}.bam".format(fq_1[:-6])
    bam_2 = "{}.bam".format(fq_2)
    if os.path.exists(Path_1):
        print("Note: The bam merging has already been performed for \
        these samples ")
        return(Path_1)
    elif os.path.exists(Path_2):
        print("Note: The bam merging has already been performed for \
        these samples ")
        return(Path_2)    
    elif not os.path.exists(fq_1):
        print("Note: The sample.bam files do not exist, the merging \
        of the files can therefore not be performed")
        return()
    else:
        cmd_bam_merge = "samtools merge -uf {} {} {}".format(Path_1, 
        bam_1, bam_2)
        res_merge = subprocess.check_output(cmd_bam_merge, shell=True)
        return(Path_1)       
    
    


if __name__ == "__main__":
    
    # Define user input:
    User_input = get_arguments()
    Samples = User_input.fastq
    Reference = User_input.reference    
    
    # Make the bowtie2 index
    bowtie2_index(Reference)
    
    # Mapping, coverging to bam, sorting, indexing and counting:
    Count_table_l = []
    Merged_bam_l = []
    
    for fq in Samples:
        # Paired samples
        if User_input.paired:
            Sam_file = bowtie2_map(Reference, fq, User_input.method)
            Bam_file = make_bam(Sam_file)
            fq_pair = find_pair(User_input.paired, fq)
            if os.path.exists(fq_pair + ".bam"):
                Merged_bam = merge_pair(fq, fq_pair)
                Merged_bam_l.append(Merged_bam)
                Sorted_bam_file = sort_bam(Merged_bam)
                index_bam(Sorted_bam_file)
                Count_table_l.append(count_bam(Sorted_bam_file))
            else:
                continue
        # Unpaired samples:
        else: 
            Sam_file = bowtie2_map(Reference, fq, User_input.method)
            Bam_file = make_bam(Sam_file)
            Sorted_bam_file = sort_bam(Bam_file)
            index_bam(Sorted_bam_file)
            Count_table_l.append(count_bam(Sorted_bam_file))

    # Parsing the count_table strings list and computing sample TPM
    TPM_l = []
    for i in range(len(Count_table_l)):
        lengths_l, counts_l = [], []
        for length, ptwy, counts in parse_count(Count_table_l[i]):
            lengths_l.append(int(length))
            counts_l.append(int(counts))
        rate_sum = calculate_rate_sum(lengths_l, counts_l)
        TPM_sample = calculate_TPM(lengths_l, counts_l, rate_sum)
        TPM_l.append(TPM_sample)

    """
    Note: 
    The formula used to calculate the TPM is:
    rate = read_count/read_length
    TPM = rate/sum(rate) * 10^6
    """

    # writing the outfile table:
    Outfile = str(User_input.output)
    with open(Outfile, "w") as o:
        o.write("...\t")
        for length, ptwy, counts in parse_count(Count_table_l[0]):
            o.write("{}\t".format(ptwy))
        o.write("\n")
        for i in range(len(TPM_l)):
            if User_input.paired:
                o.write("{}\t".format(Merged_bam_l[i][:-4]))
                for count in TPM_l[i]:
                    o.write("{}\t".format(count))
                o.write("\n")
            else:
                o.write("{}\t".format(Samples[i][:-6]))
                for count in TPM_l[i]:
                    o.write("{}\t".format(count))
                o.write("\n")
            
   
