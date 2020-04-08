#!/usr/bin/env python

"""
Name: consensus.py
Author: Xiaoyu Yu
Date: 07 April 2020

Description: Split the fasta file into multiple fasta files based on criteria set by user. For example, if the metadata file contains field lineage, the user can split the main fasta file into individual fasta files with all sequences of that lineage. The program will then create an alignment of that group of which a consensus will be built upon. The output file will consist of a single consensus file with all consensus sequences for each group divided by the trait. If a clade file is given (only work for lineages), the grouping will be collapse to the closest lineage (i.e. 1.1.2.1 collapsed to 1.1.2). The Log file will flag all sequences with no trait value and sequences that does not have a match between fasta and metadata files.

Usage:
python3 consensus.py 
--in-fasta fasta.fasta
--in-metadata metadata.csv
--index-field sequence_id
--index-column trait
--clade-file clade_file.txt
--out-fasta output_fasta

Commands: 
--in-fasta: Fasta file with sequences that needs to be splitted according to criteria to create a consensus set by user according to metadata file. (Required)
--in-metadata: Matching metadata file with same naming convention as fasta file. Contains all sequence metadata that the user wants to split the fasta file by for consensus to be created. Metadata file must be in .csv format (Required)
--index-field: The matching criteria the fasta file needs to be splitted by. (Required)
--index-column: The column with matching sequence IDs with fasta file (Default: header). (Optional)
--clade-file: Specific lineages the user wants to split by. All sub-lineages will be collapsed to the closest lineage (e.g. 1.1.2 to 1.1)
--out-folder: Output fasta file with consensus sequences of all groups based on trait (Default: consensus.fasta). (Optional)

Requirements:
BioPython (Cock et al. 2009) - http://biopython.org/wiki/Main_Page
"""

#Import Dependencies
import os
import csv
import re
import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo

#Consensus function
def create_consensus(in_fasta,in_metadata,index_field,index_column,clade_file,out_fasta):
    #Creating dictionaries for metadata, trait, sequence and consensus
    metadata_dic = {}
    trait_dic = {}
    seq_dic = {}
    consensus_dic = {}
    output_folder = os.path.dirname(out_fasta)
    log_file = open(out_fasta+".log","w")

    #Importing metadata csv to dictionary
    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    #Check if trait is in column name
    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
            sys.exit()
        else:
            metadata_dic[items[index_column]] = items[index_field.lower()]

    #Import sequence to dictionary
    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    #Check if clades are specified based on clade_file
    if os.path.isfile(clade_file):
        with open(clade_file,"r") as f:
            for line in f:
                trait_dic[line.rstrip()] = []
    else:
        for trait in metadata_dic.values():
            if trait not in trait_dic.keys():
                trait_dic[trait] = []

    #Order the traits in descending order
    trait_order = list(trait_dic.keys())
    trait_order.sort(key=lambda x: re.sub("[^A-Z0-9]", "",x),reverse=True)

    #Check if there are matching sequences to metadata names
    if len(set(metadata_dic.keys())&set(seq_dic.keys())) == 0:
        print("No matching sequence name with metadata name. Program Exit")
        sys.exit()

    #Creating traits dictionary with sequence ID and sequences
    for cluster in trait_order:
        for seq_id,phylotype in metadata_dic.items():
            cluster_type = cluster.split(".")
            cluster_length = len(cluster_type)
            phylo_type = phylotype.split(".")
            if len(phylo_type) < cluster_length:
                continue
            if phylo_type[:cluster_length] == cluster_type:
                if seq_id in seq_dic.keys():
                    trait_dic[cluster].append([seq_id,seq_dic[seq_id],phylotype])
                    del seq_dic[seq_id]

    #Check for number of sequence per trait
    for key,value in trait_dic.items():
        print(key,len(value))

    #Flag for all sequences with no matches to metadata file
    for seq in seq_dic.keys():
        log_file.write("Sequence " + seq + " did not find any matches to metadata file.\n")

    for key in trait_dic.keys():
        if len(trait_dic[key]) > 2:
            #Create fasta file for each traits with at least 2 sequences for alignment
            outfile_name = output_folder + key + ".fasta"
            outfile = open(outfile_name,"w")
            for sequences in trait_dic[key]:
                record = SeqRecord(sequences[1],id=sequences[0],description="")
                SeqIO.write(record, outfile, "fasta")
            outfile.close()
            #Create mafft alignment for each trait and removing the trait fasta file
            alignment_name = outfile_name[:-6] + "_alignment.fasta"
            align_command = "mafft " + outfile_name + " > " + alignment_name
            os.system(align_command)
            os.remove(outfile_name)
            #Creating consensus sequence using biopython to dictionary and then removing alignment fasta file
            alignment = AlignIO.read(alignment_name, 'fasta')
            consensus_name = key + "_consensus"
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus_seq = summary_align.dumb_consensus(threshold=0.0,ambiguous='N')
            consensus_dic[consensus_name] = consensus_seq
            os.remove(alignment_name)
        else:
            #Flag traits with less than two sequences
            log_file.write("Trait " + key + " does not have 2 or more sequences for an alignment to work.\n")
    log_file.close()

    #Create consensus fasta file from dictionary
    consensus_file = open(out_fasta,"w")
    for key, value in consensus_dic.items():
        record = SeqRecord(value,id=key,description="")
        SeqIO.write(record, consensus_file, "fasta")
    consensus_file.close()

#Arguments for stand alone usage
consensus = argparse.ArgumentParser(description='Collapses sequences into consensus sequences based on grouping by index column or index field')
consensus.add_argument("--in-fasta", required=True, default=None,type=str, help="<filename> [<filename> ...] one fasta files of sequences")
consensus.add_argument("--in-metadata", required=True, default=None,type=str, help="<filename> [<filename> ...] one CSV table of metadata")
consensus.add_argument("--index-field", required=True, default=None,type=str, help="<column> [<column> ...] the field(s) in the header to match the metadata")
consensus.add_argument("--index-column", required=False, default="header",type=str, help="<column> [<column> ...] the column(s) in the metadata file to use to match to the sequence")
consensus.add_argument("--clade-file", required=False, default="",type=str, help="<filename> [<filename> ...] text file including specific traits to collapse by")
consensus.add_argument("--out-fasta", required=False, default="consensus.fasta",type=str, help="<filename> output a consensus fasta file")
args = consensus.parse_args()

in_fasta = args.in_fasta
in_metadata = args.in_metadata
out_fasta  = args.out_fasta
index_field = args.index_field
index_column = args.index_column
clade_file = args.clade_file

create_consensus(in_fasta,in_metadata,index_field,index_column,clade_file,out_fasta)