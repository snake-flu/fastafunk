"""
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk)
https://github.com/cov-ert/fastafunk

This module splits sequences into fasta files based on index column or index field in metadata.

This file is part of Fastafunk. Fastafunk is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Fastafunk is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Fastafunk. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import csv
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.utils import *

def split_fasta(in_fasta,in_metadata,index_field,index_column,out_folder,log_file):
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
            metadata_dic[items[index_column]] = items[index_field.lower()]

    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

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

    for key,value in phylotype_dic.items():
        outfile = open(out_folder + key + ".fasta","w")
        for sequences in value:
            record = SeqRecord(sequences[1],id=sequences[0],description="")
            SeqIO.write(record, outfile, "fasta")
        outfile.close()
    close_handle(log_handle)