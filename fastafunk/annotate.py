"""
Name: annotate.py
Author: Rachel Colquhoun
Date: 08 April 2020
Description: Annotates the metadata file with stats calculated from the fasta file.

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


def annotate(in_fasta, in_metadata, index_column, index_field, out_fasta, out_metadata, header_delimiter, log_file):
    log_handle = get_log_handle(log_file, out_fasta)

    stats = {"length": [], "missing": [], "gaps": []}
    ids = []

    out_handle = None
    if out_fasta or not out_metadata:
        out_handle = get_out_handle(out_fasta)

    if not in_fasta:
        in_fasta = [""]

    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            record_stats = []
            for stat in stats:
                print(stat)
                result = get_stat(stat, record)
                print(result)
                record_stats.append("%s=%s" %(stat,str(result)))
                print(record_stats)
                stats[stat].append(result)
                print(stats)
            ids.append(get_index_field_from_header(record, header_delimiter, index_field))
            print(ids)
            if out_handle:
                print(record_stats)
                print(record)
                record.description += " " + " ".join(record_stats)
                print(record)
                SeqIO.write(record, out_handle, "fasta-2line")
        #close_handle(fasta_handle)

    if out_metadata:
        print("print out metadata", out_metadata)
        metadata = load_metadata(in_metadata, None, None)
        stats[index_column] = ids
        stats_data = pd.DataFrame(stats)
        metadata = add_data(stats_data, master)

        metadata_handle = get_out_handle(out_metadata)
        metadata.to_csv(out_metadata)
        #close_handle(metadata_handle)

    #close_handle(out_handle)
    #close_handle(log_handle)




