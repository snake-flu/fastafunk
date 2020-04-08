#!/usr/bin/env python

"""
Name: extract.py
Author: Xiaoyu Yu
Date: 07 April 2020

Description: Extract sequences from fasta file with matching sequence names within the metadata file. Log file will flag all sequences extracted from fasta file based on matches on metadata file.

Usage:
python3 extract.py 
--in-fasta fasta.fasta
--in-metadata metadata.csv
--out-fasta output.fasta

Commands: 
--in-fasta: Fasta file with sequences that needs to be extracted according to metadata file. (Required)
--in-metadata: Matching metadata file with same naming convention as fasta file. Contains sequences that the user wants to extract from the fasta file. Metadata file must be in .csv format (Required)
--out-fasta: Output fasta file filtered sequences extracted based on metadata file (Default: extract_by_metadata.fasta). (Optional)

Requirements:
BioPython (Cock et al. 2009) - http://biopython.org/wiki/Main_Page
"""

#Import Dependencies
import csv
from Bio import SeqIO
import argparse

#Extract function
def extract_fasta(in_fasta, in_metadata, out_fasta):

    #Creating metadata dictionary
    metadata_dictionary = {}
    with open(in_metadata) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            metadata_dictionary[row[0]] = row[1:]

    #Creating extracted fasta file
    extract_file = open(out_fasta,"w")
    log_file = open(out_fasta+".log","w")
    for record in SeqIO.parse(in_fasta, "fasta"):
        if record.id in metadata_dictionary.keys():
            SeqIO.write(record, extract_file, "fasta")
        else:
            #Flag for all extracted sequences based on matches to metadata
            log_file.write("Sequence " + record.id + " extracted due to match to metadata\n")
    extract_file.close()
    log_file.close()

#Arguments for stand alone usage
extract = argparse.ArgumentParser(description='Extracts sequences based on matches to the metadata')
extract.add_argument("--in-fasta", required=True, default=None,type=str, help="<filename> [<filename> ...] one fasta files of sequences")
extract.add_argument("--in-metadata", required=True, default=None,type=str, help="<filename> [<filename> ...] one CSV table of metadata")
extract.add_argument("--out-fasta", required=False, default="extract_by_metadata.fasta",type=str, help="<filename> output a fasta file")
args = extract.parse_args()

in_fasta = args.in_fasta
in_metadata = args.in_metadata
out_fasta  = args.out_fasta

extract_fasta(in_fasta, in_metadata, out_fasta)