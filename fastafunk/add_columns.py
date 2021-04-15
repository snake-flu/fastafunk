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

from fastafunk.metadata_reader import *
from fastafunk.utils import *

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

    all_column_names = []

    new_column_dict = {}

    new_data = MetadataReader(in_data, where_columns=where_column, index=join_on, omit_labelled_rows=False)
    if new_columns is not None and len(new_columns) > 0:
        new_data.columns = new_columns.copy()
        new_data.columns.append(join_on)
    else:
        new_columns = new_data.columns
        new_columns.remove(join_on)

    for row in new_data.reader:
        if row[join_on] in new_column_dict.keys():
            log_handle.write("Sequence " + row[join_on] + " had a duplicate in in-data and only first kept\n")
        else:
            new_column_dict[row[join_on]] = new_data.clean_row(row)
    new_data.close()

    metadata = MetadataReader(in_metadata, index=index_column, omit_labelled_rows=False)
    metadata.add_columns(new_columns)
    out_metadata_handle = open(out_metadata,"w",newline='')
    metadata.to_csv(out_metadata_handle, include_omitted=True, new_data_dict=new_column_dict)
    out_metadata_handle.close()
    metadata.close()

    log_handle.close()
