"""
Name: add_column.py
Author: Ben Jackson
Date: 15 April 2020
Description: Rewrite a metadata file with added columns(s)

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import warnings
import os
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.utils import *
from fastafunk.stats import *

def replace_with_where_columns(existing_columns, where_columns, log_handle):
    if where_columns:
        for pair in where_columns:
            column,regex = pair.split("=")
            if column in existing_columns:
                log_handle.write("Column %s already exists in in-data, ignoring %s\n" %(column, pair))
                continue
            regex = re.compile(regex)
            for existing_column in existing_columns:
                match = re.search(regex, existing_column)
                if match:
                    existing_columns[existing_columns.index(existing_column)] = column
                    log_handle.write("Renamed column %s as column %s in in-data\n" % (existing_column, column))
                    break
    return existing_columns

def add_columns(in_metadata, in_data, index_column, join_on, new_columns, out_metadata, where_column, log_file):
    """
    in_metadata - a list of metadata files to update
    in_data - a file with info used to populate new metadata columns
    index_column - column in the metadata file used to look up
                   values between the metadata and data
    join_on - column in the data file used to look up
                   values between the metadata and data
    new_columns - column in in_data that will be added to new metadata files
                  and whose values will populate the rows in new metadata, based
                  on a match between index_column, if not provided all columns added
    """
    log_handle = get_log_handle(log_file, out_fasta=False)

    join_on = join_on.lower()
    index_column = index_column.lower()
    all_column_names = []

    new_column_dict = {}
    with open(in_data, "r") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = replace_with_where_columns(reader.fieldnames, where_column, log_handle)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        if not new_columns or len(new_columns) == 0:
            new_columns = [r for r in reader.fieldnames if r!=join_on]
        data = [r for r in reader]
        new_columns = [c.lower() for c in new_columns]
    for sequence in data:
        if join_on not in sequence.keys():
            log_handle.write("Join on column not in in-data. Please re-enter a new one. Program exiting.\n")
            sys.exit()
        else:
            taxon_name = sequence[join_on]
        if taxon_name not in new_column_dict.keys():
            new_column_dict[taxon_name] = clean_dict(sequence, new_columns)
        else:
            log_handle.write("Sequence " + taxon_name + " had a duplicate in in-data and only first kept\n")

    rows = []
    null_dict = {}
    with open(in_metadata, "r") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        all_column_names = reader.fieldnames
        new_columns = [c for c in new_columns if c not in all_column_names]
        for c in new_columns:
            null_dict[c] = ''
        metadata = [clean_dict(r) for r in reader]
    for sequence in metadata:
        if index_column not in sequence.keys():
            log_handle.write("Index column not in metadata. Please re-enter a new one. Program exiting.\n")
            sys.exit()
        else:
            taxon_name = sequence[index_column]
        if taxon_name in new_column_dict.keys():
            sequence.update(new_column_dict[taxon_name])
        else:
            sequence.update(null_dict)
        rows.append(sequence)
    out_metadata_handle = open(out_metadata,"w",newline='')

    if 'unnamed: 0' in all_column_names:
        all_column_names.remove('unnamed: 0')
    if '' in all_column_names:
        all_column_names.remove('')
    all_column_names.extend(new_columns)
    f = csv.DictWriter(out_metadata_handle, fieldnames=all_column_names,lineterminator='\n')
    f.writeheader()
    f.writerows(rows)
    out_metadata_handle.close()
    log_handle.close()
