"""
Name: get_tips.py
Author: Rachel Colquhoun
Date: 15 July 2021
Description: Extract tip names corresponding to metadata.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""
import csv
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

def get_tips(in_metadata, in_tree, out_tips, low_memory):
    """
    Extract list of tips corresponding to metadata

    :param in_metadata: Matching metadata file with same naming convention as fasta file. Contains sequences that the
    user wants to extract from the fasta file. Metadata file must be in .csv format (Required)
    :param in_tree: Tree file
    :param out_tips: A tip name per line, both with and without quotes

    :return:
    """
    if in_metadata:
        metadata_dictionary = metadata_to_dict(in_metadata)
        metadata_set = set([key.lower() for key in metadata_dictionary])
        del metadata_dictionary
    else:
        sys.exit("Please specify metadata with --in-metadata")

    if in_tree:
        if low_memory:
            tree_taxon_set = wrangle_tip_labels(in_tree)
        else:
            tree_taxon_set = set(trees_to_taxa(in_tree))
        tree_taxon_set = [t.replace('"','').replace("'","") for t in tree_taxon_set]
    else:
        sys.exit("Please specify tree with --in-tree")


    out_handle = get_out_handle(out_tips)

    for tip in tree_taxon_set:
        if tip.lower() not in metadata_set:
            out_handle.write("'%s'\n" % tip)
            out_handle.write("%s\n" % tip)

    close_handle(out_handle)
