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

from fastafunk.utils import *

def subsample_fasta(in_fasta,in_metadata,index_field,index_column,group_column,out_fasta,out_metadata,log_file,
                    sample_size,target_file,select_by_max_column,select_by_min_column,exclude_uk):
    log_handle = get_log_handle(log_file, out_fasta)

    metadata = load_metadata(in_metadata, None, None)
    subsampled_metadata = subsample_metadata(metadata, group_column, sample_size, target_file, select_by_max_column,
                                             select_by_min_column, exclude_uk)

    metadata_id_key = index_column

    if not in_fasta:
        in_fasta = [""]

    out_handle = get_out_handle(out_fasta)

    to_keep = [fix_header_string(s) for s in subsampled_metadata[metadata_id_key].values]
    log_handle.write("\n#Pruned ids:\n")
    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            id_string = record.id.split('|')[0]
            if id_string in to_keep:
                SeqIO.write(record, out_handle, "fasta-2line")
                to_keep.remove(id_string)
            else:
                log_handle.write("%s\n" %record.id)
        close_handle(fasta_handle)

    if out_metadata:
        metadata_handle = get_out_handle(out_metadata)
        subsampled_metadata.to_csv(out_metadata)
        close_handle(metadata_handle)
        
    close_handle(log_handle)
    close_handle(out_handle)




