"""
Name: new.py
Author: Rachel Colquhoun
Date: 14 April 2020
Description: This module identifies the new or updated sequences in a fasta and associated metadata file.

Requires the old metadata file and a new/updated metadata file. Searches for entries in the new/updated
metadata file which have occured since the last date in the old metadata file. This subset of the metadata
is output along with associated fasta sequences. Note that it does not deduplicate, but instead outputs all
fasta sequences in input which have a header matching the subset metadata.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import csv
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.utils import *

def new_fasta(in_fasta,in_metadata,index_column,date_column,out_fasta,out_metadata,log_file,header_delimiter):
    log_handle = get_log_handle(log_file, out_fasta)

    metadata = load_new_metadata(in_metadata, date_column)
    metadata, index_column_values = get_index_column_values(metadata, index_column, header_delimiter)

    if not in_fasta:
        in_fasta = [""]

    out_handle = get_out_handle(out_fasta)

    to_keep = index_column_values
    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            id_string = record.id
            if id_string in to_keep:
                SeqIO.write(record, out_handle, "fasta-2line")
        close_handle(fasta_handle)

    if out_metadata:
        metadata_handle = get_out_handle(out_metadata)
        metadata.to_csv(out_metadata)
        close_handle(metadata_handle)
        
    close_handle(log_handle)
    close_handle(out_handle)




