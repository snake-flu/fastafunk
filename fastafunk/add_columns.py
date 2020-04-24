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
import pandas as pd

from fastafunk.utils import *
from fastafunk.stats import *


def add_columns(in_metadata, in_data, index_column, join_on, new_columns, out_metadata, log_file):
    """
    in_metadata - a list of metadata files to update
    in_data - a file with info used to populate new metadata columns
    index_column - column in the metadata file used to look up
                   values between the metadata and data
    join_on - column in the data file used to look up
                   values between the metadata and data
    new_columns - column in in_data that will be added to new metadata files
                  and whose values will populate the rows in new metadata, based
                  on a match between index_column
    """
    log_handle = get_log_handle(log_file, out_fasta=False)

    join_on = join_on.lower()
    index_column = index_column.lower()
    new_columns = [c.lower() for c in new_columns]
    all_column_names = []

    new_column_dict = {}
    with open(in_data, "r") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        data = [r for r in reader]
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
    for c in new_columns:
        null_dict[c] = ''
    with open(in_metadata, "r") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        all_column_names = reader.fieldnames
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
    f = csv.DictWriter(out_metadata_handle, fieldnames=all_column_names)
    f.writeheader()
    f.writerows(rows)
    out_metadata_handle.close()
    log_handle.close()
