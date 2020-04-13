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
    --lineage: Allow user to specify specific lineages to split by. The lineage list does not 
    need to consist of all the lineage present. All sub-lineages will collapse to the closes 
    lineage. For example --lineage A, B, B.1 will collapse all B.1* to B.1 and others to B while
    all A* will be grouped into A 

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

def split_fasta(in_fasta,in_metadata,index_field,index_column,lineage,out_folder,log_file):
    """
    Split the fasta file into multiple fasta files based on criteria set by user

    :param in_fasta: Fasta file with sequences that needs to be splitted according to criteria set by user according to
    metadata file. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Contains all sequence metadata
    that the user wants to split the fasta file by. Metadata file must be in .csv format (Required)
    :param index_field: The matching criteria the fasta file needs to be splitted by. (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
    :param lineage: Only apply to lineage specific fields. The file will consist of all the specific lineages the user wants to split by. All sub-lineages will be collapsed to the closest
    lineage (e.g. 1.1.2 to 1.1). (Optional)
    :param out_folder: Output folder for all fasta files splitted based on matching criteria (Default: ./). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    log_handle = get_log_handle(log_file, out_folder)

    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
            sys.exit()
        else:
            if items[index_column] in metadata_dic.keys():
                print("Duplicate sequences with name: " + items[index_column] + " in metadata file.", file=log_handle)
            metadata_dic[items[index_column]] = items[index_field.lower()]

    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    if len(set(metadata_dic.keys())&set(seq_dic.keys())) == 0:
        sys.exit("No matching sequence name with metadata name. Program Exit")

    if lineage != "":
        for clades in lineage:
            phylotype_dic[clades] = []

        trait_order = list(phylotype_dic.keys())
        trait_order.sort(key=lambda x: re.sub("[^A-Z0-9]", "",x),reverse=True)

        for cluster in trait_order:
            for seq_id,phylotype in metadata_dic.items():
                cluster_type = cluster.split(".")
                cluster_length = len(cluster_type)
                phylo_type = phylotype.split(".")
                if len(phylo_type) < cluster_length:
                    continue
                if phylo_type[:cluster_length] == cluster_type:
                    if seq_id in seq_dic.keys():
                        phylotype_dic[cluster].append([seq_id,seq_dic[seq_id],phylotype])
                        del seq_dic[seq_id]
    else:
        for seq,trait in metadata_dic.items():
            if trait == "":
                print("Sequence " + seq + " have an empty " + trait + " value.", file=log_handle)
            if seq not in seq_dic.keys():
                print("Sequence " + seq + " does not match metadata sequence name.", file=log_handle)
                continue
            if trait not in phylotype_dic.keys():
                phylotype_dic[trait] = []
                phylotype_dic[trait].append([seq,seq_dic[seq]])
            else:
                phylotype_dic[trait].append([seq,seq_dic[seq]])

    sort_key = sorted(phylotype_dic.keys())
    sum_values = 0
    for key in sort_key:
        print("Trait:" + key + "\t\tTotal Number:" + str(len(phylotype_dic[key])), file=log_handle)
        sum_values += len(phylotype_dic[key])
        outfile = open(out_folder + key + ".fasta","w")
        for sequences in phylotype_dic[key]:
            record = SeqRecord(sequences[1],id=sequences[0],description="")
            SeqIO.write(record, outfile, "fasta-2line")
        outfile.close()
    print("Total number of Sequences: " + str(sum_values),file=log_handle)
    close_handle(log_handle)