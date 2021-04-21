"""
Name: filter_column.py
Author: Rachel Colquhoun
Date: 24 March 2021
Description: Rewrite a metadata file with fewer rows based on filter column

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import sys


from fastafunk.utils import *


def filter_column(in_metadata, column, out_metadata, is_true, is_false, log_file):
    """
    in_metadata - a metadata file to update
    columns - columns to filter on
    out_metadata - output file
    is_true - filter out rows where column values are true
    is_false - filter out rows where column values are false
    """
    log_handle = get_log_handle(log_file, out_fasta=False)

    sep = ','
    if in_metadata.endswith('tsv'):
        sep = '\t'

    with open(in_metadata, 'r') as csv_in, open(out_metadata,'w') as csv_out:
        reader = csv.DictReader(csv_in, delimiter=sep, quotechar='\"', dialect="unix")
        writer = csv.DictWriter(csv_out, fieldnames=reader.fieldnames, delimiter=",", quotechar='\"',
                                quoting=csv.QUOTE_MINIMAL, dialect="unix")
        writer.writeheader()

        if column not in reader.fieldnames:
            sys.exit("Filter column %s not in metadata" %column)

        for row in reader:
            val = row[column].lower()
            if is_false and val in [False, "false", "f", "no", "n"]:
                continue
            if is_true and val in [True, "true", "t", "yes", "y"]:
                continue
            writer.writerow(row)

    log_handle.close()
