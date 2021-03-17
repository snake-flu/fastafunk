"""
Name: split.py
Author: Xiaoyu Yu
Date: 13 April 2020
Description: Split the fasta file into multiple fasta files based on criteria set by user.

For example, if the metadata file contains field country, the user can split the main fasta
file into individual fasta files with all sequences of that country. Log file will flag all
sequences with no trait value and sequences that does not have a match between fasta and
metadata files.

Options:
    --lineages: Allow user to specify specific lineages to split by. The lineages list does not
    need to consist of all the lineages present. All sub-lineages will collapse to the closes
    lineage. For example --lineages A, B, B.1 will collapse all B.1* to B.1 and others to B while
    all A* will be grouped into A

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import csv
import sys
import re
import json
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.utils import *

def seq_is_outgroup(seq_id, lineage_dic):
    if seq_id in lineage_dic.keys():
        return True
    return False

def get_parent(phylotype, lineages):
    parent = phylotype
    while parent != "":
        parent = ".".join(parent.split('.')[:-1])
        if parent in lineages:
            return parent
    return None

def expand_alias(phylotype, alias_dict, log_handle):
    if not phylotype or phylotype == "":
        return phylotype

    if phylotype[0] in alias_dict.keys():
        if len(phylotype) > 1:
            phylotype = alias_dict[phylotype[0]] + phylotype[1:]
        else:
            phylotype = alias_dict[phylotype[0]]
    if phylotype[0] not in ["A","B"]:
        sys.exit("Phylotype %s has no alias provided. Please update --aliases JSON" %phylotype[0])
    return phylotype

def get_clade(phylotype, ordered_lineages, alias_dict, log_handle):
    phylotype = expand_alias(phylotype, alias_dict, log_handle)
    phylo_type = phylotype.split(".")
    for cluster in ordered_lineages:
        cluster_type = cluster.split(".")
        cluster_length = len(cluster_type)
        # if the length of the phylotype is shorter than the parent, it
        # can't possibly be a match, so move on to next sequence
        if len(phylo_type) < cluster_length:
            continue
        if phylo_type[:cluster_length] == cluster_type:
            return cluster
    return "A"


def split_fasta(in_fasta,in_metadata,index_field,index_column,lineages,lineage_csv,aliases,out_prefix,log_file):
    """
    Split the fasta file into multiple fasta files based on criteria set by user

    :param in_fasta: Fasta file with sequences that needs to be splitted according to criteria set by user according to
    metadata file. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Contains all sequence metadata
    that the user wants to split the fasta file by. Metadata file must be in .csv format (Required)
    :param index_field: The matching criteria the fasta file needs to be splitted by. (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
    :param lineages: Only apply to lineages specific fields. The file will consist of all the specific lineages the user
    wants to split by. All sub-lineages will be collapsed to the closest lineages (e.g. 1.1.2 to 1.1). (Optional)
    :param lineage_csv: Only apply to lineages specific fields. The file contains two columns, named 'lineage' and
    'outgroup'. Sub-lineages are collapsed as when lineages are provided, and additionally the child outgroups are
    included in the parent file. (Optional)
    :param out_prefix: Output folder for all fasta files splitted based on matching criteria (Default: ./). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    log_handle = get_log_handle(log_file, out_prefix)

    lineage_dic = {}
    if lineage_csv:
        with open(lineage_csv) as csv_handle:
            csv_reader = csv.DictReader(csv_handle)
            for row in csv_reader:
                lineage_dic[row['outgroup']] = row['lineage']
        lineages = [lineage_dic[outgroup] for outgroup in lineage_dic.keys()]
        print("Found lineages", lineages)
    lineages.sort(key=lambda x: re.sub("[^A-Z0-9]", "", x), reverse=True)
    phylotype_counts = {}
    for lineage in lineages:
        phylotype_counts[lineage] = 0

    alias_dict = {}
    if aliases:
        print("Found alias file", aliases)
        with open(aliases, "r") as read_file:
            alias_dict = json.load(read_file)
        print("Found aliases for", alias_dict.keys())

    metadata_dic = {}
    with open(in_metadata,"r") as f:
        reader = csv.DictReader(f)
        if index_field not in reader.fieldnames:
            print("Index field name not in metadata file, please re-check metadata file.")
            sys.exit()
        if index_column not in reader.fieldnames:
            print("Index column name not in metadata file, please re-check metadata file.")
            sys.exit()
        for row in reader:
            if row[index_column] in metadata_dic:
                print("Duplicate sequences with name: " + items[index_column] + " in metadata file.", file=log_handle)
            else:
                metadata_dic[row[index_column]] = row[index_field]

    output_files = {}
    for lineage in lineages:
        filename = out_prefix + "_" + lineage + ".fasta"
        handle = open(filename, "w")
        print(filename)
        output_files[lineage] = handle

    record_dict = SeqIO.index(in_fasta, "fasta")

    for record_id in record_dict:
        if record_id not in metadata_dic:
            print("Sequence " + record_id + " does not match metadata sequence name.", file=log_handle)
            continue

        phylotype = metadata_dic[record_id]
        if phylotype in ["",None,"None"]:
            print("Sequence " + record_id + " has no lineage value.", file=log_handle)
            continue

        clade = get_clade(phylotype, lineages, alias_dict, log_handle)
        print("Add seq", record_id, "with lineage", phylotype, "to clade", clade)
        output_files[clade].write(">%s\n%s\n" % (record_id, str(record_dict[record_id].seq)))
        phylotype_counts[clade] += 1

        if seq_is_outgroup(record_id, lineage_dic):
            parent = get_parent(phylotype, lineages)
            print("Seq", record_id, "is outgroup with lineage", phylotype, "and parent lineage", str(parent))
            if parent is not None:
                output_files[parent].write(">%s\n%s\n" % (record_id, str(record_dict[record_id].seq)))
                phylotype_counts[parent] += 1

    for file in output_files:
        if hasattr(file, 'close'):
            file.close()

    sum_values = 0
    print(phylotype_counts)
    for key in phylotype_counts:
        print("Trait:" + key + "\t\tTotal Number:" + str(phylotype_counts[key]), file=log_handle)
        sum_values += phylotype_counts[key]
    print("Total number of Sequences: " + str(sum_values),file=log_handle)
    close_handle(log_handle)
