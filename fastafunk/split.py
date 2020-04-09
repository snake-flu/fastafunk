"""
Name: split.py
Author: Xiaoyu Yu
Date: 07 April 2020
Description: Split the fasta file into multiple fasta files based on criteria set by user.

For example, if the metadata file contains field country, the user can split the main fasta
file into individual fasta files with all sequences of that country. Log file will flag all
sequences with no trait value and sequences that does not have a match between fasta and
metadata files.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import csv
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.utils import *

def split_fasta(in_fasta,in_metadata,index_field,index_column,out_folder,log_file):
    """
    Split the fasta file into multiple fasta files based on criteria set by user

    :param in_fasta: Fasta file with sequences that needs to be splitted according to criteria set by user according to
    metadata file. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Contains all sequence metadata
    that the user wants to split the fasta file by. Metadata file must be in .csv format (Required)
    :param index_field: The matching criteria the fasta file needs to be splitted by. (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
    :param out_folder: Output folder for all fasta files splitted based on matching criteria (Default: ./). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    log_handle = get_log_handle(log_file, out_folder)

    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
            sys.exit()
        else:
            metadata_dic[items[index_column]] = items[index_field.lower()]

    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    for seq,trait in metadata_dic.items():
        if trait == "":
            print("Sequence " + seq + " have an empty " + trait + " value.", file=log_handle)
        if seq not in seq_dic.keys():
            print("Sequence " + seq + " does not match metadata sequence name.", file=log_handle)
            continue
        if trait not in phylotype_dic.keys():
            phylotype_dic[trait] = []
            phylotype_dic[trait].append([seq,seq_dic[seq]])
        else:
            phylotype_dic[trait].append([seq,seq_dic[seq]])

    for key,value in phylotype_dic.items():
        outfile = open(out_folder + key + ".fasta","w")
        for sequences in value:
            record = SeqRecord(sequences[1],id=sequences[0],description="")
            SeqIO.write(record, outfile, "fasta-2line")
        outfile.close()
    close_handle(log_handle)