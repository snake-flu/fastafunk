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

def merge_fasta(in_fasta, in_metadata, out_fasta, log_file):
    """
    Merges two or more fasta files avoiding duplicates based on matches to metadata

    :param in_fasta: List of fasta files with spaces in between. At least two fasta files must be inserted here. Only
    fasta files are taken as input. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Those that does not match or
    have duplicates will be flagged within the log file for post-processing. Metadata file must be in csv format
    (Required)
    :param out_fasta: Output fasta file with merged sequences from multiple inputs (Default: merged.fasta). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)
    :return:
    """
    if not in_fasta:
        in_fasta = [""]
    metadata_dictionary = metadata_to_dict(in_metadata)
    sequence_dictionary = {}
    out_handle = get_out_handle(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)

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