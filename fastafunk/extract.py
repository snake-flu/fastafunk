"""
Name: extract.py
Author: Xiaoyu Yu
Date: 07 April 2020
Description: Extract sequences from fasta file with matching sequence names within the metadata file.

Log file will flag all sequences extracted from fasta file based on matches on metadata file.


Date: 26 May 2020

Add the ability to extract sequences from a fasta file based on matched to a treefile

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import csv
from Bio import SeqIO

from fastafunk.utils import *

def extract_fasta(in_fasta, in_metadata, in_tree, out_fasta, reject_fasta, log_file):
    """
    Extract sequences from fasta file with matching sequence names within the metadata file

    :param in_fasta: Fasta file with sequences that needs to be extracted according to metadata file. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Contains sequences that the
    user wants to extract from the fasta file. Metadata file must be in .csv format (Required)
    :param out_fasta: Output fasta file filtered sequences extracted based on metadata file (Default:
    extract_by_metadata.fasta). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    if not in_fasta:
        in_fasta = [""]

    if in_metadata:
        metadata_dictionary = metadata_to_dict(in_metadata)
    else:
        metadata_dictionary = {}

    if in_tree:
        tree_taxon_set = set(trees_to_taxa(in_tree))
    else:
        tree_taxon_set = set()

    if not in_tree and not in_metadata:
        sys.exit("please specify one or both out of --in-metadata, --in-tree")

    out_handle = get_out_handle(out_fasta)
    if reject_fasta:
        reject_handle = get_out_handle(reject_fasta)
    log_handle = get_log_handle(log_file, out_fasta)

    for fasta_file in in_fasta:
        record_dict = SeqIO.index(fasta_file, "fasta")
        for record in record_dict:
            if record in metadata_dictionary.keys() or record in tree_taxon_set:
                SeqIO.write(record_dict[record], out_handle, "fasta-2line")
            elif reject_fasta:
                SeqIO.write(record_dict[record], reject_handle, "fasta-2line")
            else:
                print("Sequence " + record.id + " removed due to no match to metadata", file=log_handle)

    close_handle(reject_handle)
    close_handle(out_handle)
    close_handle(log_handle)
