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


def add_columns(in_metadata, in_data, index_column, join_on, new_columns, out_metadata, out_logfile):
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
    # log_handle = get_log_handle(log_file, out_fasta = False)

    metadata = []
    First = True
    with open(in_metadata, 'r') as f:
        for line in f:
            l = line.strip().split(',')
            if First:
                metadata_fields = l
                if not index_column in metadata_fields:
                    sys.exit(index_column + ' does not match any header in ' +  in_metadata)
                First = False
                continue
            d = {x:y for x,y in zip(metadata_fields, l)}
            metadata.append(d)

    other_data = {}
    First = True
    with open(in_data, 'r', encoding='utf-8-sig') as f:
        for line in f:
            l = line.strip().split(',')
            if First:
                header = l
                if not join_on in header:
                    sys.exit('--join-on does not match any header in the --in-data file')
                if not all([x in header for x in new_columns]):
                    sys.exit('a new column name does not match any header in the --in-data file')
                First = False
                continue
            d = {x:y for x,y in zip(header, l)}
            key = d[join_on]
            if key in other_data:
                warnings.warn(key + ' is a duplicate in ' + in_data + '[,' + join_on +'] and will be ignored')
                continue
            else:
                other_data[key] = d

    newfields = metadata_fields + new_columns

    with open(out_metadata, 'w') as f:
        f.write(','.join(newfields) + '\n')

        for dictionary in metadata:
            # this is the value in that row for index_column
            lookup = dictionary[index_column]
            # might not be an entry
            if len(lookup) == 0:
                for key in new_columns:
                    dictionary[key] = ''

                f.write(','.join([dictionary[x] for x in newfields]) + '\n')
                continue

            elif lookup in other_data:
                for key in new_columns:
                    dictionary[key] = other_data[lookup][key]

                f.write(','.join([dictionary[x] for x in newfields]) + '\n')
                continue

            else:
                for key in new_columns:
                    dictionary[key] = ''

                f.write(','.join([dictionary[x] for x in newfields]) + '\n')

    pass
