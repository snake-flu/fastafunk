"""
Name: add_column.py
Author: Ben Jackson
Date: 15 April 2020
Description: Rewrite a metadata file with added columns(s)

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.utils import *
from fastafunk.stats import *


def add_columns(in_metadata, in_data, index_column, new_columns, out_metadata, out_logfile):
    """
    in_metadata - a list of metadata files to update
    in_data - a file with info used to populate new metadata columns
    index_column - column in both metadata files and in_data file use to look up
                   values between the two
    new_columns - column in in_data that will be added to new metadata files
                  and whose values will populate the rows in new metadata, based
                  on a match between index_column
    """
    # log_handle = get_log_handle(log_file, out_fasta = False)

    metadata = load_metadata(in_metadata, None, None)

    other_data = load_dataframe(in_data, None, None)

    if not index_column in metadata.columns:
        sys.exit('index column name does not match any header in the metadata file')

    if not index_column in other_data.columns:
        sys.exit('index column name does not match any header in the in_data file')

    if not all([x in other_data.columns for x in new_columns]):
        sys.exit('a new column name does not match any header in the in_data file')

    for i in range(len(metadata)):
        # for every row in metadata:
        row = metadata.iloc[i,]
        # this is the value in that row for index_column
        lookup = row[index_column]
        # then look it up in the index_column of other_data (e.g. lineages file)
        matching_row = other_data.loc[other_data[index_column] == lookup]
        if len(matching_row) > 1:
            sys.exit('multiple entries with the same value in: ' + in_data)
        new_values = matching_row[new_columns].squeeze()

        # populate the metadata with new valus
        for key in new_values.index:
            value = new_values.at[key]
            metadata.at[i, key] = value

        metadata.to_csv(out_metadata, index = False, sep = ',')

    pass
