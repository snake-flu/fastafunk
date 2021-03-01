"""
Name: subsample.py
Author: Rachel Colquhoun
Date: 08 April 2020
Description: This module subsamples a fasta based on counts in groupings defined by
a set of index fields in metadata.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import csv
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.metadata import *
from fastafunk.utils import *

def subsample_fasta(in_fasta,in_metadata,index_field,index_column,group_column,where_field, out_fasta,out_metadata,
                    log_file,sample_size,target_file,select_by_max_column,select_by_min_column,exclude_uk):
    log_handle = get_log_handle(log_file, out_fasta)

    metadata = load_metadata(in_metadata, None, None, index_column)
    subsampled_index_column_values = subsample_metadata(metadata, group_column, sample_size, target_file, select_by_max_column,
                                             select_by_min_column, exclude_uk)

    if not in_fasta:
        in_fasta = [""]

    out_handle = get_out_handle(out_fasta)

    log_handle.write("\n#Pruned ids:\n")
    for fasta_file in in_fasta:
        #fasta_handle = get_in_handle(fasta_file)
        records = SeqIO.index(fasta_file, "fasta")
        for id_string in records:
            if id_string != "" and id_string in subsampled_index_column_values:
                SeqIO.write(records[id_string], out_handle, "fasta-2line")
            else:
                log_handle.write("%s\n" %id_string)
        #close_handle(fasta_handle)

    if out_metadata:
        metadata_handle = get_out_handle(out_metadata)
        metadata.to_csv(out_metadata, index=False)
        close_handle(metadata_handle)
        
    close_handle(log_handle)
    close_handle(out_handle)




