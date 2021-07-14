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
import re

from fastafunk.utils import *

def wrangle_tip_labels(in_tree):
    tree_taxon_set = set()
    for tree_file in in_tree:
        with open(tree_file, 'r') as tree:
            line = tree.readline()
            if "#NEXUS" in line:
                tax_labels = False
                while (not tax_labels and line):
                    line = tree.readline()
                    if "TAXLABELS" in line.upper():
                        tax_labels = True
                while (tax_labels and line):
                    if ";" in line:
                        break
                    else:
                        tree_taxon_set.add(line.strip())
                        line = tree.readline()
            else:
                while line:
                    tips = re.split(r"[,\(\);]+", line.strip())
                    tips = [t.split(":")[0] for t in tips if t.split(":")[0] != '']
                    tree_taxon_set.update(tips)
                    line = tree.readline()
        return tree_taxon_set

def extract_fasta(in_fasta, in_metadata, in_tree, out_fasta, reject_fasta, low_memory, log_file):
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
        metadata_set = set([key.lower() for key in metadata_dictionary])
        del metadata_dictionary
    else:
        metadata_set = {}

    if in_tree:
        if low_memory:
            tree_taxon_set = wrangle_tip_labels(in_tree)
        else:
            tree_taxon_set = set(trees_to_taxa(in_tree))
        tree_taxon_set = [t.lower() for t in tree_taxon_set]
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
            lowercase_record = record.lower()
            if lowercase_record in metadata_set or lowercase_record in tree_taxon_set:
                SeqIO.write(record_dict[record], out_handle, "fasta-2line")
            elif reject_fasta:
                SeqIO.write(record_dict[record], reject_handle, "fasta-2line")
            else:
                print("Sequence " + record + " removed due to no match to metadata", file=log_handle)

    if reject_fasta:
        close_handle(reject_handle)
    close_handle(out_handle)
    close_handle(log_handle)
