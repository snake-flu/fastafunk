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

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv

from fastafunk.utils import *

def merge_fasta(in_fasta, in_metadata, index_column, out_fasta, log_file):
    """
    Merges two or more fasta files avoiding duplicates based on matches to metadata

    :param in_fasta: List of fasta files with spaces in between. At least two fasta files must be inserted here. Only
    fasta files are taken as input. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Those that does not match or
    have duplicates will be flagged within the log file for post-processing. Metadata file must be in csv format
    (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
    :param out_fasta: Output fasta file with merged sequences from multiple inputs (Default: merged.fasta). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)
    :return:
    """
    if not in_fasta:
        in_fasta = [""]

    out_handle = get_out_handle(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)
    sequence_dictionary = {}
    metadata_dictionary = {}
        
    for metadata_file in in_metadata:
        with open(metadata_file,"r",encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            reader.fieldnames = [name.lower() for name in reader.fieldnames]
            metadata = [r for r in reader]

        if index_column.lower() not in reader.fieldnames:
            print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
            sys.exit()

        for items in metadata:
            if items[index_column] in metadata_dictionary.keys():
                print("Duplicate sequences with name: " + items[index_column] + " in metadata file.", file=log_handle)
            else:
                metadata_dictionary[items[index_column]] = items

    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id in metadata_dictionary.keys() and record.id not in sequence_dictionary.keys():
                sequence_dictionary[record.id] = record.seq
            elif record.id not in metadata_dictionary.keys():
                log_handle.write(record.id + " sequence is not in metadata file or the name is wrong (in file " + fasta_file + ")\n")
            elif record.id in sequence_dictionary.keys():
                log_handle.write(record.id + " is a duplicate (in file " + fasta_file + ")\n")
        close_handle(fasta_handle)

    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, out_handle, "fasta-2line")

    close_handle(out_handle)
    close_handle(log_handle)