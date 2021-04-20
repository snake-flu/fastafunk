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

from fastafunk.metadata import *
from fastafunk.utils import *
from fastafunk.stats import *


def annotate(in_fasta, in_metadata, index_column, index_field, out_fasta, out_metadata, header_delimiter,
             add_cov_id, log_file, low_memory=False):
    log_handle = get_log_handle(log_file, out_fasta)


    metadata_keys = None
    if in_metadata:
        metadata = load_metadata(in_metadata, None, None, index_column)
        metadata_keys = metadata.get_index_column_values()

    stats = {"length": [], "missing": [], "gaps": []}
    ids = []
    cov_ids = []

    out_handle = None
    if out_fasta or not out_metadata:
        out_handle = get_out_handle(out_fasta)

    if not in_fasta:
        in_fasta = [""]

    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        if low_memory:
            record_dict = SeqIO.index(fasta_handle, "fasta")
        else:
            record_dict = SeqIO.parse(fasta_handle, "fasta")
        for item in record_dict:
            if type(item) == SeqRecord:
                record = item
            else:
                record = record_dict[item]
            id = get_index_field_from_header(record, header_delimiter, index_field).split()[0]
            if metadata_keys is not None and id not in metadata_keys:
                log_handle.write("Could not find sequence header id %s in index column %s" %(id, index_column))
                continue
            record_stats = []
            for stat in stats:
                try:
                    result = get_stat(stat, record)
                except:
                    log_handle.write("Check record %s" % record.id)
                    result = 0
                record_stats.append("%s=%s" %(stat,str(result)))
                stats[stat].append(result)
            ids.append(id)
            if add_cov_id:
                cov_id = get_cov_id(record)
                cov_ids.append(cov_id)
                record_stats.append("id=%s" % cov_id)
            if out_handle:
                record.description += " " + " ".join(record_stats)
                SeqIO.write(record, out_handle, "fasta-2line")
        close_handle(fasta_handle)
    close_handle(out_handle)


    if out_metadata:
        if index_column is None or index_column == "":
            index_column = "sequence_name"
        stats[index_column] = ids
        if add_cov_id:
            stats["cov_id"] = cov_ids
        stats_metadata = Metadata(metadata_dict=stats, index=index_column)
        stats_metadata.to_csv(log_handle)
        if in_metadata:
            prev_length = len(metadata.rows)
            metadata.add_data(stats_metadata)
            new_length = len(metadata.rows)
            if new_length != prev_length:
                log_handle.write("Warning: input metadata table had length %i and now has length %i with stats"
                                 %(prev_length, new_length))
                sys.exit()
        else:
            metadata = stats_metadata
        metadata_handle = get_out_handle(out_metadata)
        metadata.to_csv(metadata_handle)
        close_handle(metadata_handle)

    close_handle(log_handle)




