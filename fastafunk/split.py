import os
import csv
import sys
import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def split_fasta(in_fasta,in_metadata,index_field,index_column,out_folder):
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    log_file = open(out_folder + "split_fasta.log","w")

    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
            sys.exit()
        else:
            metadata_dic[items[index_column]] = items[index_field.lower()]

    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    for seq,trait in metadata_dic.items():
        if trait == "":
            log_file.write("Sequence " + seq + " have an empty " + trait + " value.\n")
        if seq not in seq_dic.keys():
            log_file.write("Sequence " + seq + " does not match metadata sequence name.\n")
            continue
        if trait not in phylotype_dic.keys():
            phylotype_dic[trait] = []
            phylotype_dic[trait].append([seq,seq_dic[seq]])
        else:
            phylotype_dic[trait].append([seq,seq_dic[seq]])

    for key,value in phylotype_dic.items():
        outfile = open(out_folder + key + ".fasta","w")
        for sequences in value:
            record = SeqRecord(sequences[1],id=sequences[0],description="")
            SeqIO.write(record, outfile, "fasta")
        outfile.close()
        
    log_file.close()