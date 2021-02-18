"""
Name: fetch.py
Author: Rachel Colquhoun
Date: 18 April 2020
Description: Fetches fasta entries with a corresponding entry in a metadata file, avoiding duplicates.

Later metadata entries overwrite earlier ones.
Takes the last sequence appearance if there are duplicate entries in the in_fasta.
Only those sequences matching metadata and without a flag in an omit column will be processed into output fasta file.

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

def fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta, out_metadata, log_file, header_delimiter, low_memory):
    """
    Fetches fasta entries with a corresponding entry in a metadata file

    :param in_fasta: List of fasta files with spaces in between. At least two fasta files must be inserted here. Only
    fasta files are taken as input. (Required)
    :param in_metadata: list of matching metadata file with same naming convention as fasta file (index-column). (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: sequence_name). (Optional)
    :param out_metadata: Output metadata file with merged columns from multiple inputs (Default: None). (Optional)
    :param out_fasta: Output fasta file with merged sequences from multiple inputs (Default: stdout). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)
    :return:
    """
    log_handle = get_log_handle(log_file, out_fasta)

    metadata = load_metadata(in_metadata, filter_column, where_column)
    metadata, full_index_column_values = get_index_column_values(metadata, index_column, header_delimiter)
    subsampled_metadata = filter_by_omit_columns(metadata)
    subsampled_metadata, index_column_values = get_index_column_values(subsampled_metadata, index_column,
                                                                       header_delimiter)

    if not in_fasta:
        in_fasta = [""]

    out_handle = get_out_handle(out_fasta)
    sequence_list = []
    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        if low_memory:
            record_dict = SeqIO.index(fasta_handle, "fasta")
        else:
            record_dict = SeqIO.parse(fasta_handle, "fasta")
        for record in record_dict:
            print(type(record))
            if type(record) == SeqRecord:
                id_string = record.id
            else:
                id_string = record
            if id_string is not None and id_string in index_column_values:
                if id_string in sequence_list:
                    log_handle.write("%s is a duplicate record, keeping earliest\n" % id_string)
                elif type(record) == SeqRecord:
                    SeqIO.write(record, out_handle, "fasta-2line")
                else:
                    SeqIO.write(record_dict[id_string], out_handle, "fasta-2line")
            elif id_string is not None and id_string in full_index_column_values:
                log_handle.write("%s was marked to omit\n" %id_string)
            else:
                log_handle.write("%s has no corresponding entry in metadata table\n" %id_string)
        close_handle(fasta_handle)
        close_handle(out_handle)

    if out_metadata:
        if restrict:
            metadata = restrict_dataframe(metadata, index_column, sequence_list)
            subsampled_metadata = restrict_dataframe(subsampled_metadata, index_column, sequence_list)
        metadata_handle = get_out_handle(out_metadata)
        metadata.to_csv(out_metadata, index=False)
        close_handle(metadata_handle)

    close_handle(log_handle)
