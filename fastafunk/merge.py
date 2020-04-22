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
from fastafunk.utils import *

def clean_dict(d):
    to_delete = []
    for key in d.keys():
        if key == '':
            to_delete.append(key)
        elif "unnamed" in key:
            to_delete.append(key)
    for key in to_delete:
        del d[key]
    return d

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
    out_metadata_handle = open(out_metadata,"w",newline='')
    log_handle = get_log_handle(log_file, out_fasta)
    sequence_dictionary = {}
    metadata_dictionary = {}
    additional_rows = []
    index = index_column.lower()

    for metadata_file in in_metadata:
        if os.path.exists(metadata_file):
            with open(metadata_file,"r") as f:
                reader = csv.DictReader(f)
                reader.fieldnames = [name.lower() for name in reader.fieldnames]
                metadata = [clean_dict(r) for r in reader]
            for sequence in metadata:
                if index not in sequence.keys():
                    print("Index column not in metadata. Please re-enter a new one. Program exiting.")
                    sys.exit()
                else:
                    taxon_name = sequence[index]
                if taxon_name not in metadata_dictionary.keys():
                    metadata_dictionary[taxon_name] = sequence
                else:
                    additional_rows.append(sequence)
                    log_handle.write("Sequence " + taxon_name + " had a duplicate in metadata both kept\n")
        else:
            print("File does not exist, program exiting.")
            sys.exit()

    sequence_list = list(metadata_dictionary.keys())
    out_list = list(metadata_dictionary.values())
    all_keys = set()
    for i in out_list:
        all_keys.update(i.keys())
    all_keys.discard('unnamed: 0')
    all_keys.discard('')
    f = csv.DictWriter(out_metadata_handle, fieldnames=all_keys)
    f.writeheader()
    f.writerows(out_list)
    f.writerows(additional_rows)
    out_metadata_handle.close()
    #print(metadata_dictionary,sequence_list,out_list)

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

    if len(sequence_dictionary.keys()) == 0:
        print("There is no matching sequences to metadata. Program exiting")
        sys.exit()

    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, out_handle, "fasta-2line")

    close_handle(out_handle)
    close_handle(log_handle)
