#!/usr/bin/env python

"""
Name: split.py
Author: Xiaoyu Yu
Date: 07 April 2020

Description: Split the fasta file into multiple fasta files based on criteria set by user. For example, if the metadata file contains field country, the user can split the main fasta file into individual fasta files with all sequences of that country. Log file will flag all sequences with no trait value and sequences that does not have a match between fasta and metadata files.

Usage:
python3 split.py 
--in-fasta fasta.fasta
--in-metadata metadata.csv
--index-field trait
--index-column sequence_id
--out-fasta output_folder

Commands: 
--in-fasta: Fasta file with sequences that needs to be splitted according to criteria set by user according to metadata file. (Required)
--in-metadata: Matching metadata file with same naming convention as fasta file. Contains all sequence metadata that the user wants to split the fasta file by. Metadata file must be in .csv format (Required)
--index-field: The matching criteria the fasta file needs to be splitted by. (Required)
--index-column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
--out-folder: Output folder for all fasta files splitted based on matching criteria (Default: ./). (Optional)

Requirements:
BioPython (Cock et al. 2009) - http://biopython.org/wiki/Main_Page
"""

#Import Dependencies
import os
import csv
import sys
import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#Split function
def split_fasta(in_fasta,in_metadata,index_field,index_column,out_folder):
    #Creating dictionaries for metadata, traits and sequences
    metadata_dic = {}
    trait_dic = {}
    seq_dic = {}
    log_file = open(out_folder + "split_fasta.log","w")

    #Import metadata csv into dictionary
    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    #Check for trait existence in column names
    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
            sys.exit()
        else:
            metadata_dic[items[index_column]] = items[index_field.lower()]

    #Import sequences into sequence dictionary
    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    #Creating trait dictionary with all sequence ID and sequences
    for seq,trait in metadata_dic.items():
        if trait == "":
            #Flag for empty trait values
            log_file.write("Sequence " + seq + " have an empty " + trait + " value.\n")
        if seq not in seq_dic.keys():
            #Flag for no matches between fasta and metadata file
            log_file.write("Sequence " + seq + " does not match metadata sequence name.\n")
            continue
        if trait not in trait_dic.keys():
            trait_dic[trait] = []
            trait_dic[trait].append([seq,seq_dic[seq]])
        else:
            trait_dic[trait].append([seq,seq_dic[seq]])

    #Creating individual fasta files based on trait dictionary
    for key,value in trait_dic.items():
        outfile = open(out_folder + key + ".fasta","w")
        for sequences in value:
            record = SeqRecord(sequences[1],id=sequences[0],description="")
            SeqIO.write(record, outfile, "fasta")
        outfile.close()
        
    log_file.close()

#Arguments for stand alone usage
split = argparse.ArgumentParser(description='Splits sequences into fasta files based on index column or index field')
split.add_argument("--in-fasta", required=True, default=None,type=str, help="<filename> [<filename> ...] one fasta files of sequences")
split.add_argument("--in-metadata", required=True, default=None,type=str, help="<filename> [<filename> ...] one CSV table of metadata")
split.add_argument("--index-field", required=True, default=None,type=str, help="the field(s) in the header to match the metadata")
split.add_argument("--index-column", required=False, default="header",type=str, help="<column> [<column> ...] the column(s) in the metadata file to use to match to the sequence")
split.add_argument("--out-folder", required=False, default="./",type=str, help="<filename> output folder for fasta files split by trait")
args = split.parse_args()

in_fasta = args.in_fasta
in_metadata = args.in_metadata
out_folder  = args.out_folder
index_field = args.index_field
index_column = args.index_column

split_fasta(in_fasta,in_metadata,index_field,index_column,out_folder)