"""
Name: merge.py
Author: Xiaoyu Yu
Date: 07 April 2020
Description: Merges two or more fasta files avoiding duplicates based on matches to metadata.

Takes the first appearance according to the sequence of files within the input (--in-fasta command).
At least two fasta files must be within the input command and only those sequences matching metadata
will be processed into output fasta file.
If repeats in metadata, only first is kept

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from functools import reduce
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import sys
import os
from fastafunk.metadata_reader import *
from fastafunk.utils import *


def merge_fasta(in_fasta, in_metadata, index_column, out_metadata, out_fasta, log_file, low_memory=False):
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
        in_fasta = []

    out_handle = get_out_handle(out_fasta)
    out_metadata_handle = open(out_metadata,"w",newline='')
    log_handle = get_log_handle(log_file, out_fasta)

    sequence_list = []
    index_column_values = []
    duplicates = []
    metadata_columns = []

    for metadata_file in in_metadata:
        if os.path.exists(metadata_file):
            metadata = MetadataReader(metadata_file, None, None, index_column,omit_labelled_rows=False)
            duplicates.extend([r for r in metadata.rows if r in index_column_values])
            index_column_values.extend(metadata.rows)
            metadata_columns.extend([c for c in metadata.columns if c not in metadata_columns])
            metadata.close()
        else:
            print("File does not exist, program exiting.")
            sys.exit()

    first = True
    for metadata_file in in_metadata:
        if os.path.exists(metadata_file):
            metadata = MetadataReader(metadata_file, None, filter_columns=metadata_columns, index=index_column, omit_labelled_rows=False)
            if first:
                metadata.to_csv(out_metadata_handle, include_omitted=True, header=True)
                first = False
            else:
                metadata.omit_rows = duplicates
                metadata.to_csv(out_metadata_handle, include_omitted=False, header=False)
            metadata.close()
        else:
            print("File does not exist, program exiting.")
            sys.exit()
    for r in duplicates:
        log_handle.write("%s is a duplicate metadata record, keeping earliest\n" % r)
    out_metadata_handle.close()

    for fasta_file in in_fasta:
        if fasta_file == '':
            continue
        print("'%s'" %fasta_file)
        if low_memory:
            record_dict = SeqIO.index(fasta_file, "fasta")
        else:
            record_dict = SeqIO.parse(fasta_file, "fasta")

        for record in record_dict:
            if type(record) == SeqRecord:
                id_string = record.id
            else:
                id_string = record

            if id_string is not None:
                if not low_memory and id_string in sequence_list:
                    log_handle.write("%s is a duplicate record, keeping earliest\n" % id_string)
                elif id_string not in index_column_values:
                    log_handle.write("%s has no corresponding entry in metadata table\n" %id_string)
                elif type(record) == SeqRecord:
                    SeqIO.write(record, out_handle, "fasta-2line")
                    sequence_list.append(id_string)
                    index_column_values.remove(id_string)
                else:
                    SeqIO.write(record_dict[id_string], out_handle, "fasta-2line")
                    sequence_list.append(id_string)
                    index_column_values.remove(id_string)

    if len(sequence_list) == 0 and len(in_fasta) > 0:
        print("There is no matching sequences to metadata. Program exiting")
        sys.exit()

    close_handle(out_handle)
    close_handle(log_handle)
