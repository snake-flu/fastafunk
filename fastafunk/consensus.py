"""
Name: consensus.py
Author: Xiaoyu Yu
Date: 13 April 2020
Description: Collapses fasta sequences into consensus sequences based on criteria set by user. 

For example, if the metadata file contains field lineage, the user can split the main fasta file 
into individual fasta files with all sequences of that lineage. The program will then create an 
alignment of that group of which a consensus will be built upon. The output file will consist of 
a single consensus file with all consensus sequences for each group divided by the trait. The Log 
file will flag all sequences with no trait value and sequences that does not have a match between 
fasta and metadata files.

Options:
    --lineage: Allow user to specify specific lineages to split by. The lineage list does not need 
    to consist of all the lineage present. All sub-lineages will collapse to the closes lineage. 
    For example --lineage A, B, B.1 will collapse all B.1* to B.1 and others to B while
    all A* will be grouped into A 

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import csv
import re
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo

from fastafunk.utils import *

def create_consensus(in_fasta,in_metadata,index_field,index_column,lineage,out_fasta,log_file):
    """
    Collapses sequences into consensus sequences based on grouping by index column or index field

    :param in_fasta: Fasta file with sequences that needs to be splitted according to criteria to create a consensus
    set by user according to metadata file. (Required)
    :param in_metadata: Matching metadata file with same naming convention as fasta file. Contains all sequence
    metadata that the user wants to split the fasta file by for consensus to be created. Metadata file must be in .csv
    format (Required)
    :param index_field: The matching criteria the fasta file needs to be splitted by. (Required)
    :param index_column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
    :param lineage: Specific lineages the user wants to split by. All sub-lineages will be collapsed to the closest
    lineage (e.g. 1.1.2 to 1.1). (Optional)
    :param out_fasta: Output fasta file with consensus sequences of all groups based on trait (Default: consensus.fasta). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    consensus_dic = {}
    output_folder = os.path.dirname(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)
    log_handle.write("Output folder: %s" %output_folder)

    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            sys.exit("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
        else:
            metadata_dic[items[index_column]] = items[index_field.lower()]

    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    if len(set(metadata_dic.keys())&set(seq_dic.keys())) == 0:
        sys.exit("No matching sequence name with metadata name. Program Exit")

    print(lineage)
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

    for seq in seq_dic.keys():
        log_handle.write("Sequence " + seq + " did not find any matches to metadata file.\n")

    for key in phylotype_dic.keys():
        if len(phylotype_dic[key]) > 2:
            print("Trait:" + key + "\t\tTotal Number:" + str(len(phylotype_dic[key])), file=log_handle)
            outfile_name = output_folder + key + ".fasta"
            outfile = open(outfile_name,"w")
            for sequences in phylotype_dic[key]:
                record = SeqRecord(sequences[1],id=sequences[0],description="")
                SeqIO.write(record, outfile, "fasta-2line")
            outfile.close()
            alignment_name = outfile_name[:-6] + "_alignment.fasta"
            align_command = "mafft " + outfile_name + " > " + alignment_name
            os.system(align_command)
            os.remove(outfile_name)
            alignment = AlignIO.read(alignment_name, 'fasta')
            consensus_name = key + "_consensus"
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus_seq = summary_align.dumb_consensus(threshold=0.0,ambiguous='N')
            consensus_dic[consensus_name] = consensus_seq
            os.remove(alignment_name)
        else:
            log_handle.write("Phylotype " + key + " does not have 2 or more sequences for an alignment to work.\n")
    close_handle(log_handle)

    out_handle = get_out_handle(out_fasta)
    for key, value in consensus_dic.items():
        record = SeqRecord(value,id=key,description="")
        SeqIO.write(record, out_handle, "fasta-2line")
    close_handle(out_handle)