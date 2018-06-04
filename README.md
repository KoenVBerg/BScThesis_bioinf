# Table of contents:
- Short description
- Usage requirements
- Quickstart
- Research paper

# Short description
For the analysis of bacterial metagenomic and metatranscriptomic samples more
and more tools become available. This tool is focussed on finding the
representation of primary metabolic pathways and their related homologs in
metagenomic and metatranscriptomic samples. The tools consist of the seven
following scripts: 
1. SRA_download.py
2. fastqc_script.py
3. parse_genbank.py
4. get_pathway_fromfile.py
5. alignment.py
6. bowtie_mapping.py
7. result_processing.rmd

Each script starts with a documented section which explains the function of the
script within the whole pipeline and how the script works. For more information
on how to begin implementing the pipeling skip ahead to **Quickstart**. It is
assumed that the users knows it's way around a linux terminal and knows how to
run scripts.

# Usage requirements
## Hardware:
- Computer/server cluster consisting of at least one processing unit 
 
## Software
- Python 3.x interpreter
- Linux Ubuntu (although any distro with a terminal works)
- clustalw 2.1
- distmat EMBOSS:6.6.0.0
- fastq-dump 2.3.5
- FastQC v0.11.4
- bowtie2 2.2.6
- samtools version 0.1.19-96b5f2294a.
- Rstudio

# Quickstart
## Downloading the data
Before the actual processing and analyzing of the data the data first needs to
be downloaded. Here the data is downloaded from the online NCBI repository using
the SRA toolkit. Using the SRA_download.py script you can download the data as
follows: 
```
~$ python3 SRA_download.py -h
~$ python3 SRA_download.py -i sample_accessions.txt -s study_code -f True
```
This downloads the accession numbers mentioned in sample_accessions.txt and 
extracts the fastq files from them (mentioned with -f True). The "-h" gives a
verbose description of the usage for the script. 

## Checking read quality 
The quality assesment is performed using fastQC in fastqc_script.py as follows:
```
~$ python3 fastqc_script.py -h
~$ python3 fastqc_script.py -f fastq_files.fastq
```
This checks the read quality of the inputted fastq-files and creates a summary
table which enables you to quickly check the quality of all inputted
fastq-files.

## Obtaining pathway and homolog sequences
The pathways for this research were initially provided as genbank files.
Therefore parse_genbank.py can be used to parse the file to a fasta format. The
homolog sequences are found using multigeneblast, of which the results were also
provided. The results pages are displayed as html pages and the possibility to
pick the fasta sequence of the chosen homolog was lacking. Therefore, 
get_pathway_fromfile.py was developed. It works as follows: 
```
~$ python3 get_pathway_fromfile.py -h 
~$ python3 get_pathway_fromfile.py -n result_number -t fasta_title -i
multigeneblast_output
```
This will produce a fasta file from the multigeneblast result output by simply
picking the result number and using this script. Using alignment.py the
similarity between homologs was checked by producing a distance matrix.

## Read mapping and counting TPM
To find the representation of the pathway sequences in the metagenomic and
metatranscriptomic samples bowtie_mapping.py was developed. It uses the pathway
sequences as a reference and maps the samples against this reference. Then,
using the samtools suite the counts are calculated and converted to TPM. The
script can be used as follows: 
```
~$ python3 bowtie_mapping.py -h
~$ python3 bowtie_mapping.py -r reference -i input_fastq [-f SRA_table] -o
name_output_table -m method
```
This will produce an output table with the read counts and pathway (or
reference) sequences. First a reference index will be build to which the input
fastq files are mapped. The script includes the option for paired data, which
can be activated with the "-f" flag by specifying an SRA runtable which can be
downloaded from the NCBI repository. Also one of 4 local mapping methods can be
chosen. Use the "-h" flag for more details or take a look at the script, which
includes a elaborate explanation. 

## Results Processing
The processing of the results will be performed in an R markdown script. This is
an environment in which the commands can be run interactively to produce results
directly. Although the processing of the results varies among individuals, the
script used here, results_processing.rmd, can be used for inspiration. It uses
metagenomeSeq to normalize the results and uses ggplot2 for nice visualizations. 

# Research Paperper:
The research paper can be found as the included pdf file: Insights into metabolic pathway representation from analysis of human metagenomic and metatranscriptomic samples
