"""
Name: add_column.py
Author: Rachel Colquhoun
Date: 24 March 2021
Description: Rewrite a metadata file with fewer columns(s)

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

def drop_columns(in_metadata, columns, out_metadata, log_file):
    """
    in_metadata - a metadata file to update
    drop_columns - columns to drop
    out_metadata - output file
    """
    log_handle = get_log_handle(log_file, out_fasta=False)

    sep = ','
    if in_metadata.endswith('tsv'):
        sep = '\t'

    with open(in_metadata, 'r') as in_csv, open(out_metadata,'w') as out_csv:
        reader = csv.DictReader(in_csv, delimiter=sep)
        keep_columns = [c for c in reader.fieldnames if c not in columns]

        writer = csv.DictWriter(out_csv, fieldnames=keep_columns, delimiter=",", quotechar='\"',
                                quoting=csv.QUOTE_MINIMAL, dialect="unix")
        writer.writeheader()

        for row in reader:
            clean_row = {}
            for key in row:
                if key in keep_columns:
                    clean_row[key] = row[key]
            writer.writerow(clean_row)

    log_handle.close()
