"""
Name: merge.py
Author: Xiaoyu Yu
Date: 07 April 2020
Description: Merges two or more fasta files avoiding duplicates based on matches to metadata.

Takes the first appearance according to the sequence of files within the input (--in-fasta command).
At least two fasta files must be within the input command and only those sequences matching metadata
will be processed into output fasta file.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from functools import reduce
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import sys
import os
import pandas as pd
from fastafunk.utils import *

def merge_fasta(in_fasta, in_metadata, index_column, out_metadata, out_fasta, log_file):
    """
    Merges two or more fasta files avoiding duplicates based on matches to metadata

    :param in_fasta: List of fasta files with spaces in between. At least two fasta files must be inserted here. Only
    fasta files are taken as input. (Required)
    :param in_metadata: list of matching metadata file with same naming convention as fasta file (index-column). Those 
    that does not match or have duplicates will be flagged within the log file for post-processing. Metadata file must 
    be in csv format (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
    :param out_metadata: Output metadata file with merged columns from multiple inputs (Default: stdout). (Optional)
    :param out_fasta: Output fasta file with merged sequences from multiple inputs (Default: stdout). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)
    :return:
    """
    if not in_fasta:
        in_fasta = [""]

    if not in_metadata:
        in_metadata = [""]

    out_handle = get_out_handle(out_fasta)
    out_metadata_handle = open(out_metadata,"w",newline='', encoding='utf-8-sig')
    log_handle = get_log_handle(log_file, out_fasta)
    sequence_dictionary = {}
    metadata_dictionary = {}
    df_list = []

    for metadata_file in in_metadata:
        if os.path.exists(metadata_file):
            df = pd.read_csv(metadata_file)
            df.columns = map(str.lower, df.columns)        
            df_list.append(df)
        else:
            print("File does not exist, program exiting.")
            sys.exit()

    if len(df_list) > 1:
        group_axis = df_list[0].columns.get_loc(index_column.lower())        
        metadata_df = pd.concat(df_list, axis=group_axis, ignore_index=True, sort=False)
        metadata_df.fillna('', inplace=True)
    else:
        metadata_df = df_list[0]

    metadata_list = metadata_df.to_dict('record')
    for rows in metadata_list:
        taxon_name = rows[index_column.lower()]
        if taxon_name not in metadata_dictionary.keys():
            metadata_dictionary[taxon_name] = rows
        else:
            metadata_dictionary[taxon_name].update({k:v for k,v in rows.items() if v})
            log_handle.write("Sequence " + taxon_name + " has a duplicate in metadata and new metadata value is used\n")

    sequence_list = list(metadata_dictionary.keys())
    out_list = list(metadata_dictionary.values())
    f = csv.DictWriter(out_metadata_handle, fieldnames=out_list[0].keys())
    f.writeheader()
    f.writerows(out_list)
    out_metadata_handle.close()

    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id in sequence_list and record.id not in sequence_dictionary.keys():
                sequence_dictionary[record.id] = record.seq
            elif record.id not in sequence_list:
                log_handle.write(record.id + " sequence is not in metadata file or the name is wrong (in file " + fasta_file + ")\n")
            elif record.id in sequence_dictionary.keys():
                log_handle.write(record.id + " is a duplicate (in file " + fasta_file + ")\n")
        close_handle(fasta_handle)

    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, out_handle, "fasta-2line")

    close_handle(out_handle)
    close_handle(log_handle)
