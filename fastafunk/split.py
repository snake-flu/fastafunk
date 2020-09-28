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

def seq_is_outgroup(seq_id, lineage_dic):
    if seq_id in lineage_dic.keys():
        return True
    return False

def get_parent(phylotype, lineage):
    parent = phylotype
    while parent != "":
        parent = ".".join(parent.split('.')[:-1])
        if parent in lineage:
            return parent
    return None

def split_fasta(in_fasta,in_metadata,index_field,index_column,lineage,lineage_csv,out_folder,log_file):
    """
    Split the fasta file into multiple fasta files based on criteria set by user

    :param in_fasta: Fasta file with sequences that needs to be splitted according to criteria set by user according to
    metadata file. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Contains all sequence metadata
    that the user wants to split the fasta file by. Metadata file must be in .csv format (Required)
    :param index_field: The matching criteria the fasta file needs to be splitted by. (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
    :param lineage: Only apply to lineage specific fields. The file will consist of all the specific lineages the user
    wants to split by. All sub-lineages will be collapsed to the closest lineage (e.g. 1.1.2 to 1.1). (Optional)
    :param lineage_csv: Only apply to lineage specific fields. The file contains two columns, named 'lineage' and
    'outgroup'. Sub-lineages are collapsed as when lineages are provided, and additionally the child outgroups are
    included in the parent file. (Optional)
    :param out_folder: Output folder for all fasta files splitted based on matching criteria (Default: ./). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    lineage_dic = {}
    log_handle = get_log_handle(log_file, out_folder)

    with open(in_metadata,"r") as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
        print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
        sys.exit()
    for items in metadata:
        if items[index_column] in metadata_dic.keys():
            print("Duplicate sequences with name: " + items[index_column] + " in metadata file.", file=log_handle)
        else:
            metadata_dic[items[index_column]] = items[index_field.lower()]

    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    if len(set(metadata_dic.keys())&set(seq_dic.keys())) == 0:
        sys.exit("No matching sequence name with metadata name. Program Exit")

    if lineage_csv:
        with open(lineage_csv) as csv_handle:
            csv_reader = csv.DictReader(csv_handle)
            for row in csv_reader:
                lineage_dic[row['outgroup']] = row['lineage']
        lineage = [lineage_dic[outgroup] for outgroup in lineage_dic.keys()]
        print("Found lineages", lineage)

    if lineage != "":
        # clades are the lineage rows in lineage_splits
        for clades in lineage:
            phylotype_dic[clades] = []
        # print(phylotype_dic)
        trait_order = list(phylotype_dic.keys())
        # trait order is a list of the parent lineages in lineage_splits
        # print(trait_order)

        trait_order.sort(key=lambda x: re.sub("[^A-Z0-9]", "",x),reverse=True)
        # print(trait_order)

        for cluster in trait_order:
            # cluster is parent lineage in lineage_splits
            # phylotype is A, B.1, B.1.X, etc. - pangolin assigned lineage
            for seq_id,phylotype in metadata_dic.items():

                # Another hack by Ben
                #   - if phylotype is C.X/D.X, then set phylotype to B.1.1.1.X/B.1.1.25.X
                if phylotype[0] == "C":
                    if len(phylotype) > 1:
                        phylotype = "B.1.1.1" + phylotype[1:]
                    else:
                        phylotype = "B.1.1.1"

                if phylotype[0] == "D":
                    if len(phylotype) > 1:
                        phylotype = "B.1.1.25" + phylotype[1:]
                    else:
                        phylotype = "B.1.1.25"

                # print(seq_id,phylotype)
                cluster_type = cluster.split(".")
                # print(cluster_type)
                cluster_length = len(cluster_type)
                # print(cluster_length)
                phylo_type = phylotype.split(".")
                # print(phylo_type)
                # if the length of the phylotype is shorter than the parent, it
                # can't possibly be a match, so move on to next sequence
                if len(phylo_type) < cluster_length:
                    continue
                # print(phylo_type[:cluster_length])

                # if there's a match, then do the correct stuff:
                if phylo_type[:cluster_length] == cluster_type:
                    if seq_id in seq_dic.keys():
                        phylotype_dic[cluster].append([seq_id,seq_dic[seq_id],phylotype])
                        if seq_is_outgroup(seq_id, lineage_dic):
                            parent = get_parent(phylotype, lineage)
                            print("Seq", seq_id, "is outgroup with lineage" , phylotype, "and parent lineage", str(parent))
                            if parent is not None:
                                phylotype_dic[parent].append([seq_id, seq_dic[seq_id], phylotype])
                        del seq_dic[seq_id]

                # Ben's addition which is a bit of a hack - if A is the only lineage defined
                # in trait_orders, then assign all Bs to A:
                else:
                    if len(trait_order) == 1 and trait_order == ['A']:
                        if phylo_type[:1] != ['A']:
                            if seq_id in seq_dic.keys():
                                phylotype_dic[cluster].append([seq_id,seq_dic[seq_id],phylotype])
                                if seq_is_outgroup(seq_id, lineage_dic):
                                    parent = get_parent(phylotype, lineage)
                                    print("Seq", seq_id, "is outgroup with lineage" , phylotype, "and parent lineage", str(parent))
                                    if parent is not None:
                                        phylotype_dic[parent].append([seq_id, seq_dic[seq_id], phylotype])
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
