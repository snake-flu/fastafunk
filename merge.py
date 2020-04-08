#!/usr/bin/env python

"""
Name: merge.py
Author: Xiaoyu Yu
Date: 07 April 2020

Description: Merges two or more fasta files avoiding duplicates based on matches to metadata (takes the first appearance according to the sequence of files within the input (--in-fasta command)). At least two fasta files must be within the input command and only those sequences matching metadata will be processed into output fasta file.

Usage:
python3 merge.py 
--in-fasta fasta1.fasta fasta2.fasta ... 
--in-metadata metadata.csv
--out-fasta output.fasta

Commands: 
--in-fasta: List of fasta files with spaces in between. At least two fasta files must be inserted here. Only fasta files are taken as input. (Required)
--in-metadata: Matching metadata file with same naming convention as fasta file. Those that does not match or have duplicates will be flagged within the log file for post-processing. Metadata file must be in csv format(Required)
--out-fasta: Output fasta file with merged sequences from multiple inputs (Default: merged.fasta). (Optional)

Requirements:
BioPython (Cock et al. 2009) - http://biopython.org/wiki/Main_Page
"""

#Import Dependencies
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import argparse

#Merge Function with imports from file
def merge_fasta(in_fasta, in_metadata, out_fasta):
    #Creating metadata and sequence dictionaries.
    metadata_dictionary = {}
    with open(in_metadata) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            metadata_dictionary[row[0]] = row[1:]
    
    #Creating merged sequence dictionary for output
    log_file = open(out_fasta+".log","w")    
    sequence_dictionary = {}
    for fasta_file in in_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in metadata_dictionary.keys() and record.id not in sequence_dictionary.keys():
                sequence_dictionary[record.id] = record.seq
            elif record.id not in metadata_dictionary.keys():
                #Flagging for none matching sequence names
                log_file.write(record.id + " sequence is not in metadata file or the name is wrong (in file " + fasta_file + ")\n")
            elif record.id in sequence_dictionary.keys():
                #Flagging for duplicates
                log_file.write(record.id + " is a duplicate (in file " + fasta_file + ")\n")

    #Output merged sequences to file from dictionary
    merged_file = open(out_fasta,"w")
    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, merged_file, "fasta")

    merged_file.close()
    log_file.close()

#Arguments for stand alone usage
merge = argparse.ArgumentParser(description='Merges two fasta files avoiding duplicates based on matches to metadata')
merge.add_argument("--in-fasta", nargs='+', required=True, default=None,type=str, help="<filename> [<filename> ...] two or more fasta files of sequences")
merge.add_argument("--in-metadata", required=True, default=None,type=str, help="<filename> [<filename> ...] one CSV table of metadata")
merge.add_argument("--out-fasta", required=False, default="merged.fasta",type=str, help="<filename> output a fasta file")
args = merge.parse_args()

in_fasta = args.in_fasta
in_metadata = args.in_metadata
out_fasta  = args.out_fasta

merge_fasta(in_fasta, in_metadata, out_fasta)